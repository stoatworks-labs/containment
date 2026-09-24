/**
 * Containment — browser demo.
 *
 * **The solver on this page is the plugin's solver, run for real.** Every
 * piece of GLSL in `source/Shaders.cpp` is copied below unedited, and the
 * programs are assembled from those pieces exactly as `SourceFor()` joins
 * them: the reduction, the clock, MUSCL-Hancock's predictor, the HLLD face
 * fluxes with GLM, the update with EGLM, the dual-energy switch, the floors,
 * cooling and the absorbing margin, then ignite and sources, the emission,
 * the glow, the multigrid for the field lines and the composite. The state is
 * four RGBA32F textures ping-ponged, as in the plugin, and dt lives on the
 * GPU. `demo/tools/check_shaders.py` proves every piece and every assembly
 * character for character, and `tools/verify.sh` runs it.
 *
 * What is a hand port, and is checked by nothing but a reader:
 * `Controls.cpp` (every 0..1 conversion), the CPU half of `Physics.cpp`
 * (`MakeCoils`, `VacuumPotential`, `ChooseGrid`, the PCG generator and the
 * Ornstein-Uhlenbeck `Drive`), `Presets.h` and the override in
 * `ContainmentPlugin::P`, and the frame sequence of `ProcessOpenGL`: the
 * substep plan from last frame's read-back speed plus 10% (capped at 48), the
 * synchronous plan after an Ignite or a pellet, the coils' spin and quench,
 * the carried coil change, the events, `EnsureBuffers`, `Glow`'s Gaussian
 * weights and `Potential`'s V-cycle schedule. The parameter declarations come
 * from the constructor in `Containment.cpp`.
 *
 * **What the page runs by default is not the plugin's default grid.** See
 * DEMO_DETAIL below: the Detail control defaults to 128 here, where the
 * plugin's constructor sets 256, because at 256 a browser cannot keep up with
 * ~22 substeps of seven dependent passes a frame. That is said in the banner,
 * in the disclosure and in AGENTS.md, and 256, 512 and 1024 are all still on
 * the dropdown, running the same solver.
 *
 * **Ignite and Pellet are buttons.** The plugin declares them FF_TYPE_EVENT
 * and acts on the rising edge; the kit has no event type, so each is a
 * boolean the renderer takes and releases on the frame it acts — as millpond's
 * Drop does — which is why it blinks. Quench is FF_TYPE_BOOLEAN in the plugin
 * and a plain toggle here.
 *
 * **Audio is absent.** Audio, Audio Heat and Audio Pellets read Resolume's FFT
 * buffer; a browser has none. They are left off the panel rather than shown
 * dead, and the removal is exact: with no spectrum the plugin's analyser
 * reports a level of 0 and never fires, so Audio Heat adds no heat and Audio
 * Pellets drops no pellet — which is what the page does.
 */

import { mountDemo } from './vendor/demo.js';
import { Program, PassBuffer, GLError, bindTexture } from './vendor/gl.js';

//===========================================================================
// The shaders. Copied from source/Shaders.cpp. Do not edit here.
//===========================================================================

const VERSION = `#version 410 core
`;

const VERTEX = `
layout( location = 0 ) in vec4 vPosition;
layout( location = 1 ) in vec2 vUV;
out vec2 uv;
void main()
{
	gl_Position = vPosition;
	uv          = vUV;
}
`;

const COMMON = `
uniform ivec2 GridSize;   //cells, the margin included
uniform float Dx;         //cell size, frame heights (square cells)
uniform ivec2 FrameOrigin;//the frame's first cell: the margin (Open), else 0
uniform ivec2 FrameCells; //the frame's own cells

//The clip's coordinates at a cell. The margin reads the clip's edge.
vec2 clipUV( ivec2 cell )
{
	return clamp( ( vec2( cell - FrameOrigin ) + 0.5 ) / vec2( FrameCells ), 0.0, 1.0 );
}
uniform float Gamma;      //5/3; 2 for Brio-Wu
uniform int BoundaryMode; //0 open, 1 conducting wall, 2 periodic (harness only)
uniform int UseFloors;    //1: rho and p are clamped to the floors below
uniform float RhoFloor;
uniform float PFloor;
uniform float AmbientDensity; //the background the frame sits in (Open)
uniform float AmbientPressure;

//The coils: ( x, y, I / 2 pi ) each, frame heights, the rotation and any quench
//already applied on the CPU; and the uniform guide field.
uniform int CoilCount;
uniform vec3 CoilData[ 12 ];
uniform float GuideBz;

//A state, or a flux, in primitive form, as four vec4s so one limiter serves
//all thirteen quantities:
//  q0 ( rho, u, v, w )  q1 ( p, Bx, By, Bz )  q2 ( psi, K, 0, 0 )  q3 ( r, g, b, chi )
//K = p / rho^gamma is the entropy the dual-energy switch carries (see update).
struct Q
{
	vec4 q0;
	vec4 q1;
	vec4 q2;
	vec4 q3;
};

vec2 cellCentre( ivec2 cell )
{
	return ( vec2( cell ) + 0.5 ) * Dx;
}

//= mirrored in Physics.cpp, VacuumField()
//Infinite line currents along z: B = sum I' ( -(y - yk), x - xk ) / d^2.
vec3 vacuumField( vec2 x )
{
	vec2 b = vec2( 0.0 );
	for( int k = 0; k < CoilCount; ++k )
	{
		vec2 d = x - CoilData[ k ].xy;
		b += CoilData[ k ].z * vec2( -d.y, d.x ) / dot( d, d );
	}
	return vec3( b, GuideBz );
}

//= mirrored in Physics.cpp, VacuumPotential()
float vacuumPotential( vec2 x )
{
	float a = 0.0;
	for( int k = 0; k < CoilCount; ++k )
	{
		vec2 d = x - CoilData[ k ].xy;
		a -= CoilData[ k ].z * 0.5 * log( dot( d, d ) );
	}
	return a;
}

Q toPrim( vec4 A, vec4 B, vec4 C, vec4 D )
{
	Q w;
	float r = A.x;
	if( UseFloors == 1 )
		r = max( r, RhoFloor );
	vec3 v = A.yzw / r;
	float p = ( Gamma - 1.0 ) * ( B.x - 0.5 * dot( A.yzw, v ) - 0.5 * dot( B.yzw, B.yzw ) );
	if( UseFloors == 1 )
		p = max( p, PFloor );
	w.q0 = vec4( r, v );
	w.q1 = vec4( p, B.yzw );
	w.q2 = vec4( C.x, C.z / r, 0.0, 0.0 );
	w.q3 = D / r;
	return w;
}

void toCons( Q w, out vec4 A, out vec4 B, out vec4 C, out vec4 D )
{
	float r = w.q0.x;
	vec3 v  = w.q0.yzw;
	vec3 b  = w.q1.yzw;
	A = vec4( r, r * v );
	B = vec4( w.q1.x / ( Gamma - 1.0 ) + 0.5 * r * dot( v, v ) + 0.5 * dot( b, b ), b );
	//The entropy is always made consistent with the pressure when a state is
	//built from primitives: a pass that sets p (Ignite, a heating) sets K.
	C = vec4( w.q2.x, 0.0, r * w.q1.x / pow( r, Gamma ), 0.0 );
	D = r * w.q3;
}

//The fast magnetosonic speed along a direction whose field component is bn.
//Written as ( a^2 - b^2 )^2 + 4 a^2 bt^2 under the root, which cannot go
//negative by rounding the way ( a^2 + b^2 )^2 - 4 a^2 bn^2 can.
float fastSpeed( float r, float p, vec3 B, float bn )
{
	float a2  = Gamma * p / r;
	float b2  = dot( B, B ) / r;
	float bt2 = max( dot( B, B ) - bn * bn, 0.0 ) / r;
	float d   = a2 - b2;
	return sqrt( 0.5 * ( a2 + b2 + sqrt( d * d + 4.0 * a2 * bt2 ) ) );
}

//Signal speeds of a cell: ( |u| + cfx + |v| + cfy, max of the two ). The CFL
//step and GLM's c_h both come from the sum.
vec2 signalSpeeds( Q w )
{
	float sx = abs( w.q0.y ) + fastSpeed( w.q0.x, w.q1.x, w.q1.yzw, w.q1.y );
	float sy = abs( w.q0.z ) + fastSpeed( w.q0.x, w.q1.x, w.q1.yzw, w.q1.z );
	return vec2( sx + sy, max( sx, sy ) );
}

//A state, from any four textures holding one, with the boundary applied to
//cells off the grid. Ghosts are built in primitive form.
Q fetchPrim( sampler2D TA, sampler2D TB, sampler2D TC, sampler2D TD, ivec2 cell )
{
	ivec2 n      = GridSize;
	bool inside  = cell.x >= 0 && cell.y >= 0 && cell.x < n.x && cell.y < n.y;
	ivec2 source = cell;

	//Per axis: 0 open, 1 wall, 2 periodic, 3 outflow (harness), 4 periodic
	//in x with walls in y (harness: the Rayleigh-Taylor slab).
	bool periodicX = BoundaryMode == 2 || BoundaryMode == 4;
	bool periodicY = BoundaryMode == 2;
	bool wallX     = BoundaryMode == 1;
	bool wallY     = BoundaryMode == 1 || BoundaryMode == 4;
	bool offX      = cell.x < 0 || cell.x >= n.x;
	bool offY      = cell.y < 0 || cell.y >= n.y;
	if( !inside )
	{
		if( periodicX )
			source.x = ( cell.x + n.x ) % n.x;//cell >= -2 and n >= 8, so never negative
		else if( wallX )
			source.x = cell.x < 0 ? -1 - cell.x : ( cell.x >= n.x ? 2 * n.x - 1 - cell.x : cell.x );
		else
			source.x = clamp( cell.x, 0, n.x - 1 );
		if( periodicY )
			source.y = ( cell.y + n.y ) % n.y;
		else if( wallY )
			source.y = cell.y < 0 ? -1 - cell.y : ( cell.y >= n.y ? 2 * n.y - 1 - cell.y : cell.y );
		else
			source.y = clamp( cell.y, 0, n.y - 1 );
	}

	Q w = toPrim( texelFetch( TA, source, 0 ), texelFetch( TB, source, 0 ), texelFetch( TC, source, 0 ),
	              texelFetch( TD, source, 0 ) );

	if( ( offX && wallX ) || ( offY && wallY ) )
	{
		//Reflecting: the plasma is mirrored. The field is NOT mirrored: a
		//field line crossing a wall is continuous through it, and a mirrored
		//one has a kink at the face that launches a wave every step (it did,
		//and blew up). The ghost's field is the coils' field there plus the
		//mirrored cell's departure from the coils' field: an exact
		//equilibrium when nothing has happened, and the plasma's own
		//compression carried through. No slip: the whole velocity turns
		//round, so the face is at rest. A slip wall that field lines cross
		//cannot also be a perfect conductor -- E_t = v_t B_n is not zero --
		//and the mismatch with the zero flux of B through the face made a
		//current sheet on the wall that heated and flung the plasma (|v| 80
		//in 50 frames). Line-tied and at rest, E = 0 at the face and every
		//override in the flux pass is exact.
		w.q0.yzw = -w.q0.yzw;
		w.q2.x   = -w.q2.x;
		w.q1.yzw += vacuumField( cellCentre( cell ) ) - vacuumField( cellCentre( source ) );
	}
	else if( !inside && BoundaryMode == 0 )
	{
		//Open: the frame is a window onto a bigger bottle. The plasma carries
		//on as it was (zero gradient), so whatever leaks out leaves -- but the
		//ghost never holds LESS than the ambient background. A pure
		//zero-gradient ghost let the tenuous background drain out along the
		//field lines that cross the edge, leaving beta ~ 0 holes that only
		//the floors held up (millions of clamps, then a runaway). The field
		//is the coils' own at the ghost, and GLM's psi leaves.
		w.q0.x   = max( w.q0.x, AmbientDensity );
		w.q1.x   = max( w.q1.x, AmbientPressure );
		w.q2.y   = w.q1.x / pow( w.q0.x, Gamma );
		w.q1.yzw += vacuumField( cellCentre( cell ) ) - vacuumField( cellCentre( source ) );
		w.q2.x   = 0.0;
	}
	return w;
}
`;

const RIEMANN = `
uniform int Solver; //0 HLLD (Miyoshi & Kusano 2005), 1 HLL
uniform int UseGLM; //1: Dedner's hyperbolic-parabolic divergence cleaning

float gCh = 1.0;    //GLM's wave speed this substep, set by main() from the clock

//Monotonised central: minmod( 2 dl, 2 dr, ( dl + dr ) / 2 ).
vec4 mcLimit( vec4 l, vec4 r )
{
	vec4 m = min( min( 2.0 * abs( l ), 2.0 * abs( r ) ), 0.5 * abs( l + r ) );
	return sign( l + r ) * m * vec4( greaterThan( l * r, vec4( 0.0 ) ) );
}

Q slopeOf( Q m, Q c, Q p )
{
	Q s;
	s.q0 = mcLimit( c.q0 - m.q0, p.q0 - c.q0 );
	s.q1 = mcLimit( c.q1 - m.q1, p.q1 - c.q1 );
	s.q2 = mcLimit( c.q2 - m.q2, p.q2 - c.q2 );
	s.q3 = mcLimit( c.q3 - m.q3, p.q3 - c.q3 );
	return s;
}

Q addScaled( Q a, Q b, float k )
{
	Q s;
	s.q0 = a.q0 + k * b.q0;
	s.q1 = a.q1 + k * b.q1;
	s.q2 = a.q2 + k * b.q2;
	s.q3 = a.q3 + k * b.q3;
	return s;
}

//One side of a face, rotated so n is the face normal and t, w are tangential.
struct Side
{
	float r, vn, vt, vw, p, bn, bt, bw, psi, K;
	vec4 c;
};

Side rotate( Q w, int dir )
{
	Side s;
	s.r   = w.q0.x;
	s.p   = w.q1.x;
	s.vw  = w.q0.w;
	s.bw  = w.q1.w;
	s.psi = w.q2.x;
	s.K   = w.q2.y;
	s.c   = w.q3;
	if( dir == 0 )
	{
		s.vn = w.q0.y; s.vt = w.q0.z;
		s.bn = w.q1.y; s.bt = w.q1.z;
	}
	else
	{
		s.vn = w.q0.z; s.vt = w.q0.y;
		s.bn = w.q1.z; s.bt = w.q1.y;
	}
	return s;
}

//A flux in rotated form: m ( rho, mn, mt, mw ), e ( E, Bn, Bt, Bw ).
struct Flux
{
	vec4 m;
	vec4 e;
	float psi;
	float ent;//the flux of rho K
	vec4 c;
};

float energyOf( Side s, float bn )
{
	return s.p / ( Gamma - 1.0 ) + 0.5 * s.r * ( s.vn * s.vn + s.vt * s.vt + s.vw * s.vw )
	       + 0.5 * ( bn * bn + s.bt * s.bt + s.bw * s.bw );
}

void physical( Side s, float bn, float E, out vec4 fm, out vec4 fe )
{
	float pt = s.p + 0.5 * ( bn * bn + s.bt * s.bt + s.bw * s.bw );
	float vb = s.vn * bn + s.vt * s.bt + s.vw * s.bw;
	fm = vec4( s.r * s.vn, s.r * s.vn * s.vn + pt - bn * bn, s.r * s.vn * s.vt - bn * s.bt, s.r * s.vn * s.vw - bn * s.bw );
	fe = vec4( ( E + pt ) * s.vn - bn * vb, 0.0, s.bt * s.vn - bn * s.vt, s.bw * s.vn - bn * s.vw );
}

//HLLD, following Miyoshi & Kusano (2005) section by section, with the two
//degenerate cases handled as they recommend: a tangential field that cannot
//be solved for (denominator ~ 0) keeps its side's value, and a normal field
//~ 0 collapses the Alfven fans onto the contact.
void hlld( Side L, Side R, float bn, out vec4 fm, out vec4 fe )
{
	float bn2 = bn * bn;
	float EL  = energyOf( L, bn );
	float ER  = energyOf( R, bn );
	float ptL = L.p + 0.5 * ( bn2 + L.bt * L.bt + L.bw * L.bw );
	float ptR = R.p + 0.5 * ( bn2 + R.bt * R.bt + R.bw * R.bw );

	vec3 bL = vec3( bn, L.bt, L.bw );
	vec3 bR = vec3( bn, R.bt, R.bw );
	float cfL = fastSpeed( L.r, L.p, bL, bn );
	float cfR = fastSpeed( R.r, R.p, bR, bn );
	//(67): the outermost signal speeds.
	float SL = min( L.vn, R.vn ) - max( cfL, cfR );
	float SR = max( L.vn, R.vn ) + max( cfL, cfR );

	vec4 FLm, FLe, FRm, FRe;
	physical( L, bn, EL, FLm, FLe );
	physical( R, bn, ER, FRm, FRe );

	vec4 ULm = vec4( L.r, L.r * L.vn, L.r * L.vt, L.r * L.vw );
	vec4 ULe = vec4( EL, 0.0, L.bt, L.bw );
	vec4 URm = vec4( R.r, R.r * R.vn, R.r * R.vt, R.r * R.vw );
	vec4 URe = vec4( ER, 0.0, R.bt, R.bw );

	if( SL >= 0.0 )
	{
		fm = FLm; fe = FLe;
		return;
	}
	if( SR <= 0.0 )
	{
		fm = FRm; fe = FRe;
		return;
	}
	if( Solver == 1 )
	{
		//HLL: one intermediate state, no contact, no Alfven waves.
		float inv = 1.0 / ( SR - SL );
		fm = ( SR * FLm - SL * FRm + SL * SR * ( URm - ULm ) ) * inv;
		fe = ( SR * FLe - SL * FRe + SL * SR * ( URe - ULe ) ) * inv;
		return;
	}

	float sdL = SL - L.vn;
	float sdR = SR - R.vn;
	float den = sdR * R.r - sdL * L.r;
	//(38) the contact's speed, (41) the total pressure across it.
	float SM  = ( sdR * R.r * R.vn - sdL * L.r * L.vn - ptR + ptL ) / den;
	float pts = ( sdR * R.r * ptL - sdL * L.r * ptR + L.r * R.r * sdR * sdL * ( R.vn - L.vn ) ) / den;
	float sdmL = SL - SM;
	float sdmR = SR - SM;

	//(43)-(47) the outer star states.
	float rLs = L.r * sdL / sdmL;
	float rRs = R.r * sdR / sdmR;
	float vtLs = L.vt, vwLs = L.vw, btLs = L.bt, bwLs = L.bw;
	float dL = L.r * sdL * sdmL - bn2;
	if( abs( dL ) >= 1.0e-5 * pts )
	{
		float inv = 1.0 / dL;
		float mv  = bn * ( SM - L.vn ) * inv;
		vtLs = L.vt - L.bt * mv;
		vwLs = L.vw - L.bw * mv;
		float mb = ( L.r * sdL * sdL - bn2 ) * inv;
		btLs = L.bt * mb;
		bwLs = L.bw * mb;
	}
	float vtRs = R.vt, vwRs = R.vw, btRs = R.bt, bwRs = R.bw;
	float dR = R.r * sdR * sdmR - bn2;
	if( abs( dR ) >= 1.0e-5 * pts )
	{
		float inv = 1.0 / dR;
		float mv  = bn * ( SM - R.vn ) * inv;
		vtRs = R.vt - R.bt * mv;
		vwRs = R.vw - R.bw * mv;
		float mb = ( R.r * sdR * sdR - bn2 ) * inv;
		btRs = R.bt * mb;
		bwRs = R.bw * mb;
	}
	//(48) their energies.
	float vbL  = L.vn * bn + L.vt * L.bt + L.vw * L.bw;
	float vbLs = SM * bn + vtLs * btLs + vwLs * bwLs;
	float ELs  = ( sdL * EL - ptL * L.vn + pts * SM + bn * ( vbL - vbLs ) ) / sdmL;
	float vbR  = R.vn * bn + R.vt * R.bt + R.vw * R.bw;
	float vbRs = SM * bn + vtRs * btRs + vwRs * bwRs;
	float ERs  = ( sdR * ER - ptR * R.vn + pts * SM + bn * ( vbR - vbRs ) ) / sdmR;

	vec4 ULsm = vec4( rLs, rLs * SM, rLs * vtLs, rLs * vwLs );
	vec4 ULse = vec4( ELs, 0.0, btLs, bwLs );
	vec4 URsm = vec4( rRs, rRs * SM, rRs * vtRs, rRs * vwRs );
	vec4 URse = vec4( ERs, 0.0, btRs, bwRs );

	//(51) the Alfven waves.
	float sqL = sqrt( rLs );
	float sqR = sqrt( rRs );
	float SLs = SM - abs( bn ) / sqL;
	float SRs = SM + abs( bn ) / sqR;

	if( SLs >= 0.0 )
	{
		fm = FLm + SL * ( ULsm - ULm );
		fe = FLe + SL * ( ULse - ULe );
		return;
	}
	if( SRs <= 0.0 )
	{
		fm = FRm + SR * ( URsm - URm );
		fe = FRe + SR * ( URse - URe );
		return;
	}

	if( 0.5 * bn2 < 1.0e-5 * pts )
	{
		//No normal field: the rotational discontinuities sit on the contact
		//and the double-star states do not exist.
		if( SM >= 0.0 )
		{
			fm = FLm + SL * ( ULsm - ULm );
			fe = FLe + SL * ( ULse - ULe );
		}
		else
		{
			fm = FRm + SR * ( URsm - URm );
			fe = FRe + SR * ( URse - URe );
		}
		return;
	}

	//(59)-(63) the inner states.
	float inv  = 1.0 / ( sqL + sqR );
	float sg   = sign( bn );
	float vtss = ( sqL * vtLs + sqR * vtRs + ( btRs - btLs ) * sg ) * inv;
	float vwss = ( sqL * vwLs + sqR * vwRs + ( bwRs - bwLs ) * sg ) * inv;
	float btss = ( sqL * btRs + sqR * btLs + sqL * sqR * ( vtRs - vtLs ) * sg ) * inv;
	float bwss = ( sqL * bwRs + sqR * bwLs + sqL * sqR * ( vwRs - vwLs ) * sg ) * inv;
	float vbss = SM * bn + vtss * btss + vwss * bwss;

	if( SM >= 0.0 )
	{
		float ELss = ELs - sqL * ( vbLs - vbss ) * sg;
		vec4 Um = vec4( rLs, rLs * SM, rLs * vtss, rLs * vwss );
		vec4 Ue = vec4( ELss, 0.0, btss, bwss );
		fm = FLm + SLs * Um - ( SLs - SL ) * ULsm - SL * ULm;
		fe = FLe + SLs * Ue - ( SLs - SL ) * ULse - SL * ULe;
	}
	else
	{
		float ERss = ERs + sqR * ( vbRs - vbss ) * sg;
		vec4 Um = vec4( rRs, rRs * SM, rRs * vtss, rRs * vwss );
		vec4 Ue = vec4( ERss, 0.0, btss, bwss );
		fm = FRm + SRs * Um - ( SRs - SR ) * URsm - SR * URm;
		fe = FRe + SRs * Ue - ( SRs - SR ) * URse - SR * URe;
	}
}

//The flux through a face from its two reconstructed sides. GLM first
//(Dedner et al. 2002, eq. 42): the normal field and psi at the face come from
//their own 2x2 linear Riemann problem, and HLLD then runs with that normal
//field on both sides.
Flux faceFlux( Q wl, Q wr, int dir )
{
	Side L = rotate( wl, dir );
	Side R = rotate( wr, dir );
	float bn, psim;
	if( UseGLM == 1 )
	{
		bn   = 0.5 * ( L.bn + R.bn ) - 0.5 * ( R.psi - L.psi ) / gCh;
		psim = 0.5 * ( L.psi + R.psi ) - 0.5 * gCh * ( R.bn - L.bn );
	}
	else
	{
		bn   = 0.5 * ( L.bn + R.bn );
		psim = 0.0;
	}

	Flux F;
	hlld( L, R, bn, F.m, F.e );
	F.e.y = psim;                                     //d Bn / dt + d psi / dn = 0
	F.psi = UseGLM == 1 ? gCh * gCh * bn : 0.0;       //d psi / dt + c_h^2 d Bn / dn = 0
	//The picture and the ball marker ride on the mass flux, upwind.
	F.c   = F.m.x * ( F.m.x >= 0.0 ? L.c : R.c );
	F.ent = F.m.x * ( F.m.x >= 0.0 ? L.K : R.K );
	return F;
}

//The physical flux of one state, in rotated form, GLM terms included: for
//the Hancock half step, which differences a cell's own two face states.
Flux cellFlux( Q w, int dir )
{
	Side s = rotate( w, dir );
	Flux F;
	physical( s, s.bn, energyOf( s, s.bn ), F.m, F.e );
	F.e.y = UseGLM == 1 ? s.psi : 0.0;
	F.psi = UseGLM == 1 ? gCh * gCh * s.bn : 0.0;
	F.c   = s.r * s.vn * s.c;
	F.ent = s.r * s.vn * s.K;
	return F;
}

//Back to the textures' layout: A ( rho, mx, my, mz ), B ( E, Bx, By, Bz ).
void unrotate( Flux F, int dir, out vec4 fA, out vec4 fB )
{
	if( dir == 0 )
	{
		fA = F.m;
		fB = F.e;
	}
	else
	{
		fA = vec4( F.m.x, F.m.z, F.m.y, F.m.w );
		fB = vec4( F.e.x, F.e.z, F.e.y, F.e.w );
	}
}
`;

const FORCES = `
uniform float Curvature;     //g_eff at the reference temperature, pointing away from GravityCentre
uniform float CurvatureTemp; //the reference temperature; <= 0: g_eff uniform (the harness's RT slab)
uniform vec2 GravityCentre;  //frame heights: the ball's centre
uniform float GravityCore;   //the radius inside which g_eff tapers to 0
uniform int DriveCount;
uniform vec4 DriveModes[ 12 ];//( kx, ky, a, b ): phi += a cos( k.x ) + b sin( k.x )
uniform float DriveWindow;   //the stirring's reach about the ball, frame heights

vec2 bodyForce( vec2 x, float T )
{
	vec2 d  = x - GravityCentre;
	float r = length( d );
	//The curvature drift of a real bottle's bad-curvature region goes as
	//T / R_c, so the effective gravity is proportional to the temperature:
	//g_eff = Curvature T / T_ref. Its force density is then Curvature p /
	//T_ref -- continuous across an isobaric interface, and jumping where the
	//pressure does, at the ball's edge, which is where interchange is driven.
	//A uniform g_eff (the first version) also pulled the cold background
	//outwards, stratified it by e^( g L / T_bg ) ~ 70 across the frame and
	//emptied the space round the ball. Radial, smoothed through the centre
	//where the direction is undefined: r / sqrt( r^2 + core^2 ).
	float scale = CurvatureTemp > 0.0 ? T / CurvatureTemp : 1.0;
	vec2 g = Curvature * scale * d / sqrt( r * r + GravityCore * GravityCore );

	//The stirring force, curl( phi W z ) = ( d( phi W ) / dy, -d( phi W ) / dx ):
	//divergence-free by construction, however phi and W are chosen. W is a
	//Gaussian window about the ball, so the stirring keeps the ball alive
	//without churning the tenuous background (which, stirred, opened
	//rarefied pockets whose Alfven speed set the whole frame's time step).
	if( DriveCount > 0 )
	{
		float W   = exp( -0.5 * dot( d, d ) / ( DriveWindow * DriveWindow ) );
		vec2 gradW = -W * d / ( DriveWindow * DriveWindow );
		float phi = 0.0;
		vec2 gradPhi = vec2( 0.0 );
		for( int m = 0; m < DriveCount; ++m )
		{
			vec4 mode = DriveModes[ m ];
			float ph  = dot( mode.xy, x );
			float s   = sin( ph );
			float c   = cos( ph );
			phi += mode.z * c + mode.w * s;
			gradPhi += ( -mode.z * s + mode.w * c ) * mode.xy;
		}
		vec2 grad = W * gradPhi + phi * gradW;
		g += vec2( grad.y, -grad.x );
	}
	return g;
}
`;

const REDUCE_MAIN = `
uniform sampler2D Source;
uniform ivec2 SourceSize;
uniform int FromState;//1: Source is the state's C texture ( psi, speed, rho K, flags )
uniform int Block;    //texels per side each output reduces
out vec4 fragColor;

//( max summed signal speed, floors fired, cells on the entropy branch ) of
//one texel.
vec3 reading( ivec2 p )
{
	vec4 t = texelFetch( Source, p, 0 );
	//flags = floor fired + 2 x on the entropy branch
	vec3 v = FromState == 1 ? vec3( t.y, mod( t.w, 2.0 ), floor( t.w * 0.5 ) ) : t.xyz;
	//A cell that has gone non-finite stops the clock rather than vanishing
	//from the maximum (max() with a NaN is undefined).
	if( isnan( v.x ) || isinf( v.x ) )
		v.x = 3.0e38;
	return v;
}

void main()
{
	ivec2 origin = ivec2( gl_FragCoord.xy ) * Block;
	vec3 acc = vec3( 0.0 );
	for( int j = 0; j < Block; ++j )
		for( int i = 0; i < Block; ++i )
		{
			ivec2 p = origin + ivec2( i, j );
			if( p.x >= SourceSize.x || p.y >= SourceSize.y )
				continue;
			vec3 v = reading( p );
			acc.x = max( acc.x, v.x );
			acc.y += v.y;
			acc.z += v.z;
		}
	fragColor = vec4( acc, 0.0 );
}
`;

const CLOCK_MAIN = `
uniform sampler2D Reduced;  //the last reduction level ( speed, floors, entropy cells )
uniform ivec2 ReducedSize;
uniform sampler2D PrevClock;//2x1: ( dt, done, c_h, floors ), ( entropy-branch cells, 0, 0, 0 )
uniform int Substep;
uniform int Substeps;
uniform float Target;       //simulated time this frame should cover
uniform float Courant;
uniform float Dx;
uniform float Eta;
out vec4 fragColor;

void main()
{
	vec3 r = vec3( 0.0 );
	for( int j = 0; j < ReducedSize.y; ++j )
		for( int i = 0; i < ReducedSize.x; ++i )
		{
			vec3 t = texelFetch( Reduced, ivec2( i, j ), 0 ).xyz;
			r.x = max( r.x, t.x );
			r.y += t.y;
			r.z += t.z;
		}
	vec4 prev    = texelFetch( PrevClock, ivec2( 0 ), 0 );
	//The second texel tallies the cells on the dual-energy switch, summed
	//over the frame's substeps.
	if( int( gl_FragCoord.x ) == 1 )
	{
		vec4 tally = texelFetch( PrevClock, ivec2( 1, 0 ), 0 );
		fragColor  = vec4( ( Substep == 0 ? 0.0 : tally.x ) + r.z, 0.0, 0.0, 0.0 );
		return;
	}
	float done   = Substep == 0 ? 0.0 : prev.y;
	float floors = ( Substep == 0 ? 0.0 : prev.w ) + r.y;

	float dt = Courant * Dx / max( r.x, 1.0e-20 );
	//Explicit diffusion in 2-D is stable to dx^2 / ( 4 eta ).
	if( Eta > 0.0 )
		dt = min( dt, 0.2 * Dx * Dx / Eta );
	dt = min( dt, max( Target - done, 0.0 ) / float( Substeps - Substep ) );
	//c_h, GLM's cleaning speed, is HALF the summed signal speed. Its waves
	//cross cells in both directions at once, and the unsplit step is stable
	//for the SUM of the two Courant numbers <= 1: at c_h = sum that sum is
	//2 x 0.8 and psi blows up in a few steps (it did). At sum / 2 it is 0.8,
	//and c_h is still about the fastest wave's speed wherever the two
	//directions' speeds are comparable.
	fragColor = vec4( dt, done + dt, max( 0.5 * r.x, 1.0e-6 ), floors );
}
`;

const PREDICT_MAIN = `
uniform sampler2D StateA;
uniform sampler2D StateB;
uniform sampler2D StateC;
uniform sampler2D StateD;
uniform sampler2D ClockTexture;

layout( location = 0 ) out vec4 outA;
layout( location = 1 ) out vec4 outB;
layout( location = 2 ) out vec4 outC;
layout( location = 3 ) out vec4 outD;

void main()
{
	ivec2 cell = ivec2( gl_FragCoord.xy );
	vec4 clk   = texelFetch( ClockTexture, ivec2( 0 ), 0 );
	float dt   = clk.x;
	gCh        = clk.z;

	Q c  = fetchPrim( StateA, StateB, StateC, StateD, cell );
	Q xm = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( -1, 0 ) );
	Q xp = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 1, 0 ) );
	Q ym = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 0, -1 ) );
	Q yp = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 0, 1 ) );

	Q sx = slopeOf( xm, c, xp );
	Q sy = slopeOf( ym, c, yp );

	//Hancock: the cell's own face states, a half step on by the difference of
	//their physical fluxes.
	vec4 aR, bR, aL, bL, aT, bT, aB, bB;
	Flux fR = cellFlux( addScaled( c, sx, 0.5 ), 0 );
	Flux fL = cellFlux( addScaled( c, sx, -0.5 ), 0 );
	Flux fT = cellFlux( addScaled( c, sy, 0.5 ), 1 );
	Flux fB = cellFlux( addScaled( c, sy, -0.5 ), 1 );
	unrotate( fR, 0, aR, bR );
	unrotate( fL, 0, aL, bL );
	unrotate( fT, 1, aT, bT );
	unrotate( fB, 1, aB, bB );

	float k = 0.5 * dt / Dx;
	vec4 A = texelFetch( StateA, cell, 0 );
	vec4 B = texelFetch( StateB, cell, 0 );
	vec4 C = texelFetch( StateC, cell, 0 );
	vec4 D = texelFetch( StateD, cell, 0 );
	A -= k * ( ( aR - aL ) + ( aT - aB ) );
	B -= k * ( ( bR - bL ) + ( bT - bB ) );
	C.x -= k * ( ( fR.psi - fL.psi ) + ( fT.psi - fB.psi ) );
	C.z -= k * ( ( fR.ent - fL.ent ) + ( fT.ent - fB.ent ) );
	D -= k * ( ( fR.c - fL.c ) + ( fT.c - fB.c ) );

	//The body force, a half step.
	vec2 g = bodyForce( cellCentre( cell ), c.q1.x / c.q0.x );
	A.yz += 0.5 * dt * c.q0.x * g;
	B.x  += 0.5 * dt * c.q0.x * dot( c.q0.yz, g );

	outA = A;
	outB = B;
	outC = C;
	outD = D;
}
`;

const FLUX_MAIN = `
uniform sampler2D StateA;
uniform sampler2D StateB;
uniform sampler2D StateC;
uniform sampler2D StateD;
uniform sampler2D StarA;
uniform sampler2D StarB;
uniform sampler2D StarC;
uniform sampler2D StarD;
uniform sampler2D ClockTexture;
uniform int Direction;      //0: the x faces, ( nx + 1 ) x ny; 1: the y faces, nx x ( ny + 1 )
uniform float Eta;          //resistivity

//( rho, mx, my, mz ) and ( E, Bx, By, Bz ) fluxes; ( psi, rho K, the normal
//field's GLM flux, 0 ); the picture's.
layout( location = 0 ) out vec4 outA;
layout( location = 1 ) out vec4 outB;
layout( location = 2 ) out vec4 outC;
layout( location = 3 ) out vec4 outD;

Q star( ivec2 cell )
{
	return fetchPrim( StarA, StarB, StarC, StarD, cell );
}

vec3 fieldAt( ivec2 cell )
{
	return fetchPrim( StateA, StateB, StateC, StateD, cell ).q1.yzw;
}

void main()
{
	ivec2 face = ivec2( gl_FragCoord.xy );
	ivec2 e    = Direction == 0 ? ivec2( 1, 0 ) : ivec2( 0, 1 );
	gCh        = texelFetch( ClockTexture, ivec2( 0 ), 0 ).z;

	//The face between cells L and R. Computed ONCE, here, and read by both
	//cells in the update: the flux one cell loses is bit for bit the flux
	//its neighbour gains.
	ivec2 R = face;
	ivec2 L = face - e;
	Q l2 = star( L - e );
	Q l1 = star( L );
	Q r1 = star( R );
	Q r2 = star( R + e );
	Flux F = faceFlux( addScaled( l1, slopeOf( l2, l1, r1 ), 0.5 ), addScaled( r1, slopeOf( l1, r1, r2 ), -0.5 ),
	                   Direction );

	int along   = Direction == 0 ? face.x : face.y;
	int extent  = Direction == 0 ? GridSize.x : GridSize.y;
	bool walled = BoundaryMode == 1 || ( BoundaryMode == 4 && Direction == 1 );
	bool wall   = walled && ( along == 0 || along == extent );
	if( wall )
	{
		F.m.x = 0.0;
		F.e   = vec4( 0.0 );
		F.c   = vec4( 0.0 );
		F.ent = 0.0;
	}

	vec4 fA, fB;
	unrotate( F, Direction, fA, fB );

	//Resistivity: E = eta J added to the fluxes of B and of E. In 2.5-D,
	//J = ( dBz/dy, -dBz/dx, dBy/dx - dBx/dy ), from the state at the start of
	//the step (first order in time for the diffusion, the slow part).
	if( Eta > 0.0 && !wall )
	{
		ivec2 t = ivec2( e.y, e.x );//across the face
		vec3 bL = fieldAt( L ), bR = fieldAt( R );
		vec3 bLp = fieldAt( L + t ), bLm = fieldAt( L - t );
		vec3 bRp = fieldAt( R + t ), bRm = fieldAt( R - t );
		float h  = 1.0 / Dx;
		vec3 b   = 0.5 * ( bL + bR );
		if( Direction == 0 )
		{
			//Jy = -dBz/dx, Jz = dBy/dx - dBx/dy. d By/dt = d Ez/dx, so the x
			//flux of By is -eta Jz and of Bz is eta Jy; energy, eta ( J x B )_x.
			float jy = -( bR.z - bL.z ) * h;
			float jz = ( bR.y - bL.y ) * h - 0.25 * ( ( bLp.x - bLm.x ) + ( bRp.x - bRm.x ) ) * h;
			fB += vec4( Eta * ( jy * b.z - jz * b.y ), 0.0, -Eta * jz, Eta * jy );
		}
		else
		{
			//Jx = dBz/dy, Jz = dBy/dx - dBx/dy. d Bx/dt = -d Ez/dy, so the y
			//flux of Bx is eta Jz and of Bz is -eta Jx; energy, eta ( J x B )_y.
			float jx = ( bR.z - bL.z ) * h;
			float jz = 0.25 * ( ( bLp.y - bLm.y ) + ( bRp.y - bRm.y ) ) * h - ( bR.x - bL.x ) * h;
			fB += vec4( Eta * ( jz * b.x - jx * b.z ), Eta * jz, 0.0, -Eta * jx );
		}
	}

	outA = fA;
	outB = fB;
	outC = vec4( F.psi, F.ent, F.e.y, 0.0 );
	outD = F.c;
}
`;

const UPDATE_MAIN = `
uniform sampler2D StateA;
uniform sampler2D StateB;
uniform sampler2D StateC;
uniform sampler2D StateD;
uniform sampler2D StarA;
uniform sampler2D StarB;
uniform sampler2D ClockTexture;
uniform sampler2D FluxXA;
uniform sampler2D FluxXB;
uniform sampler2D FluxXC;
uniform sampler2D FluxXD;
uniform sampler2D FluxYA;
uniform sampler2D FluxYB;
uniform sampler2D FluxYC;
uniform sampler2D FluxYD;
uniform float CoolingRate;  //Lambda
uniform float CoolingFloorT;//no cooling below this temperature
uniform float GLMAlpha;
uniform int UseEntropy;     //the dual-energy switch
uniform float EntropySwitch;//p below this fraction of kinetic + magnetic: use K
uniform float SpongeEFolds; //0: no absorbing layer; else what the fastest wave loses crossing it
uniform float SpongePower;  //the rate rises as depth^SpongePower

layout( location = 0 ) out vec4 outA;
layout( location = 1 ) out vec4 outB;
layout( location = 2 ) out vec4 outC;
layout( location = 3 ) out vec4 outD;

void main()
{
	ivec2 cell = ivec2( gl_FragCoord.xy );
	vec4 clk   = texelFetch( ClockTexture, ivec2( 0 ), 0 );
	float dt   = clk.x;
	gCh        = clk.z;

	ivec2 right = cell + ivec2( 1, 0 );
	ivec2 top   = cell + ivec2( 0, 1 );
	float k = dt / Dx;
	vec4 A = texelFetch( StateA, cell, 0 );
	vec4 B = texelFetch( StateB, cell, 0 );
	vec4 C = texelFetch( StateC, cell, 0 );
	vec4 D = texelFetch( StateD, cell, 0 );
	A -= k * ( ( texelFetch( FluxXA, right, 0 ) - texelFetch( FluxXA, cell, 0 ) )
	           + ( texelFetch( FluxYA, top, 0 ) - texelFetch( FluxYA, cell, 0 ) ) );
	B -= k * ( ( texelFetch( FluxXB, right, 0 ) - texelFetch( FluxXB, cell, 0 ) )
	           + ( texelFetch( FluxYB, top, 0 ) - texelFetch( FluxYB, cell, 0 ) ) );
	vec4 cxl = texelFetch( FluxXC, cell, 0 ), cxr = texelFetch( FluxXC, right, 0 );
	vec4 cyb = texelFetch( FluxYC, cell, 0 ), cyt = texelFetch( FluxYC, top, 0 );
	D -= k * ( ( texelFetch( FluxXD, right, 0 ) - texelFetch( FluxXD, cell, 0 ) )
	           + ( texelFetch( FluxYD, top, 0 ) - texelFetch( FluxYD, cell, 0 ) ) );

	//The cleaning does no work on the gas. GLM's psi moves the NORMAL field
	//through the faces (the x faces' Bx, the y faces' By), and in plain
	//GLM-MHD that changes B^2/2 while E stays put -- so the pressure pays
	//for it, and where beta is 1e-3 it goes negative (it did, at the open
	//boundary, every run). Dedner et al. (2002)'s EGLM energy source,
	//dE/dt = -B.grad(psi), in the form that is exact per step: the change in
	//magnetic energy the psi fluxes caused is put back into E.
	if( UseGLM == 1 )
	{
		vec2 dPsi = -k * vec2( cxr.z - cxl.z, cyt.z - cyb.z );
		vec2 bOld = B.yz - dPsi;
		B.x += 0.5 * ( dot( B.yz, B.yz ) - dot( bOld, bOld ) );
	}
	C.x -= k * ( ( cxr.x - cxl.x ) + ( cyt.x - cyb.x ) );
	C.z -= k * ( ( cxr.y - cxl.y ) + ( cyt.y - cyb.y ) );

	//The half-step state, for the sources.
	Q c = toPrim( texelFetch( StarA, cell, 0 ), texelFetch( StarB, cell, 0 ), vec4( 0.0 ), vec4( 0.0 ) );

	//The body force at the half step (midpoint rule).
	vec2 g = bodyForce( cellCentre( cell ), c.q1.x / c.q0.x );
	A.yz += dt * c.q0.x * g;
	B.x  += dt * dot( c.q0.x * c.q0.yz, g );

	//GLM's parabolic half: psi decays (Mignone & Tzeferacos 2010).
	if( UseGLM == 1 )
		C.x *= exp( -GLMAlpha * gCh * dt / Dx );

	//-------------------------------------------------------------------
	// The pressure. Where beta is tiny, p = ( gamma - 1 )( E - kinetic -
	// magnetic ) is a small difference of big numbers, and the scheme's own
	// truncation error in B^2/2 is bigger than p itself: at beta 1e-3 (the
	// background beside the coils) it went negative every step. There, the
	// pressure comes from the entropy K carried alongside (the dual-energy
	// switch of Ryu et al. 1993 / Balsara & Spicer 1999), and E is made to
	// agree. Everywhere else the energy decides and K is resynchronised to
	// it, so shocks heat as they should. Each switch is counted.
	//-------------------------------------------------------------------
	float r      = max( A.x, UseFloors == 1 ? RhoFloor : 1e-30 );
	float kin    = 0.5 * dot( A.yzw, A.yzw ) / r;
	float mag    = 0.5 * dot( B.yzw, B.yzw );
	float p      = ( Gamma - 1.0 ) * ( B.x - kin - mag );
	float onEntropy = 0.0;
	if( UseEntropy == 1 && p < EntropySwitch * ( kin + mag ) )
	{
		p         = max( C.z / r, 0.0 ) * pow( r, Gamma );
		onEntropy = 1.0;
	}

	//Bremsstrahlung, dE/dt = -Lambda rho^2 sqrt(T), integrated exactly at
	//fixed density: d sqrt(p) / dt = -( gamma - 1 ) Lambda rho^1.5 / 2. It
	//stops at CoolingFloorT, the background's temperature.
	if( CoolingRate > 0.0 )
	{
		float pMin = r * CoolingFloorT;
		if( p > pMin )
		{
			float root = sqrt( p ) - 0.5 * ( Gamma - 1.0 ) * CoolingRate * r * sqrt( r ) * dt;
			p = max( root > 0.0 ? root * root : 0.0, pMin );
		}
	}

	//The floors, counted.
	float fired = 0.0;
	if( UseFloors == 1 )
	{
		if( !( A.x >= RhoFloor ) )
		{
			float scale = A.x > 0.0 ? RhoFloor / A.x : 0.0;
			A.yzw *= scale;
			D *= scale;
			A.x   = RhoFloor;
			kin   = 0.5 * dot( A.yzw, A.yzw ) / A.x;
			fired = 1.0;
		}
		if( !( p >= PFloor ) )
		{
			p     = PFloor;
			fired = 1.0;
		}
	}
	//E and K agree with the pressure chosen, whichever way it was chosen.
	if( onEntropy > 0.0 || CoolingRate > 0.0 || fired > 0.0 )
		B.x = p / ( Gamma - 1.0 ) + kin + mag;
	C.z = A.x * p / pow( A.x, Gamma );

	//-------------------------------------------------------------------
	// Open's absorbing layer, in the margin outside the frame and nowhere
	// else. The frame is a window onto a bigger bottle, full of the ambient
	// plasma at rest in the coils' own field; in the margin the plasma is
	// relaxed towards exactly that, at a rate rising as the square of the
	// depth, and sized so that the fastest wave (c_h) loses SpongeEFolds
	// e-foldings on the way in -- and as many again on the way back. The
	// relaxation is exact for the step, so no rate is too stiff. The picture
	// is kept; the ball marker fades with the rest.
	//-------------------------------------------------------------------
	if( SpongeEFolds > 0.0 && FrameOrigin.x > 0 )
	{
		vec2 x     = cellCentre( cell );
		vec2 lo    = vec2( FrameOrigin ) * Dx;
		vec2 hi    = vec2( FrameOrigin + FrameCells ) * Dx;
		vec2 out2  = max( max( lo - x, x - hi ), vec2( 0.0 ) );
		float w    = float( FrameOrigin.x ) * Dx;
		float s    = clamp( max( out2.x, out2.y ) / w, 0.0, 1.0 );
		if( s > 0.0 )
		{
			//A wave at v crossing the layer loses rate * w / ( 3 v ).
			float rate = ( SpongePower + 1.0 ) * SpongeEFolds * gCh / w * pow( s, SpongePower );
			float f    = exp( -rate * dt );
			Q q        = toPrim( A, B, C, D );
			q.q0.x     = AmbientDensity + ( q.q0.x - AmbientDensity ) * f;
			q.q0.yzw  *= f;
			q.q1.x     = AmbientPressure + ( q.q1.x - AmbientPressure ) * f;
			vec3 bv    = vacuumField( x );
			q.q1.w     = bv.z + ( q.q1.w - bv.z ) * f;
			q.q2.x    *= f;
			q.q3.a    *= f;
			toCons( q, A, B, C, D );
		}
	}

	Q w  = toPrim( A, B, C, D );
	outA = A;
	outB = B;
	outC = vec4( C.x, signalSpeeds( w ).x, C.z, fired + 2.0 * onEntropy );
	outD = D;
}
`;

const IGNITE_MAIN = `
uniform sampler2D InputTexture;
uniform vec2 MaxUV;
uniform vec2 BallCentre;
uniform float BallRadius;
uniform int ProfileKind;  //0 Gaussian, 1 top hat
uniform float BallPressure;
uniform float ClipHeats;
uniform float BgDensity;
uniform float BgPressure;

layout( location = 0 ) out vec4 outA;
layout( location = 1 ) out vec4 outB;
layout( location = 2 ) out vec4 outC;
layout( location = 3 ) out vec4 outD;

void main()
{
	ivec2 cell = ivec2( gl_FragCoord.xy );
	vec2 x     = cellCentre( cell );
	vec4 clip  = texture( InputTexture, clipUV( cell ) * MaxUV );

	float r = length( x - BallCentre );
	float profile;
	if( ProfileKind == 0 )
		profile = exp( -( r * r ) / ( BallRadius * BallRadius ) );
	else
		//Clamped: tanh( x ) is ( e^x - e^-x ) / ( e^x + e^-x ) on this GPU,
		//inf / inf = NaN for |x| past ~44, and a top hat's far field is
		//hundreds of cells out -- every cell outside the ball came out NaN.
		profile = 0.5 - 0.5 * tanh( clamp( ( r - BallRadius ) / ( 0.75 * Dx ), -20.0, 20.0 ) );

	//Bright parts of the picture push harder. At ClipHeats 1 a white pixel
	//starts at twice the ball's pressure and a black one at none of it.
	float luma = dot( clip.rgb, vec3( 0.2126, 0.7152, 0.0722 ) );
	float heat = 1.0 - ClipHeats + 2.0 * ClipHeats * luma;

	Q w;
	w.q0 = vec4( BgDensity + ( 1.0 - BgDensity ) * profile, 0.0, 0.0, 0.0 );
	w.q1 = vec4( BgPressure + ( BallPressure * heat - BgPressure ) * profile, vacuumField( x ) );
	w.q1.x = max( w.q1.x, BgPressure * ( 1.0 - profile ) + PFloor );
	w.q2 = vec4( 0.0 );
	w.q3 = vec4( clamp( clip.rgb, 0.0, 1.0 ), profile );

	vec4 A, B, C, D;
	toCons( w, A, B, C, D );
	outA = A;
	outB = B;
	outC = vec4( C.x, signalSpeeds( w ).x, C.z, 0.0 );
	outD = D;
}
`;

const SOURCES_MAIN = `
uniform sampler2D StateA;
uniform sampler2D StateB;
uniform sampler2D StateC;
uniform sampler2D StateD;
uniform sampler2D InputTexture;
uniform vec2 MaxUV;
uniform float FeedFraction; //this frame's step of the tracer towards the clip
uniform float ClipHeatRate; //dp for a white pixel inside the ball, this frame
uniform float AudioHeatRate;//dp inside the ball, this frame
uniform float FuelFraction; //this frame's step of the footprint towards its ignition profile
uniform vec2 BallCentre;
uniform float BallRadius;
uniform int ProfileKind;
uniform float BallPressure;
uniform float BgDensity;
uniform float BgPressure;
uniform int PelletCount;    //pellets landing this frame
uniform vec2 PelletCentre;
uniform float PelletDensity;
uniform float PelletRadius;
//The coils as they were last frame. With a conducting vessel (Wall) their
//slow changes -- a turn, a quench -- are carried straight through the vessel
//and the plasma: B += B_vac( now ) - B_vac( then ), which is curl-free and
//divergence-free, so it adds no current and no monopole, and the plasma has
//to follow the field it finds itself in. (Driving the change through the
//wall's faces instead needs E_z = -dA/dt at a no-slip, line-tied wall, where
//ideal MHD cannot carry it: it made a current sheet and a runaway.)
uniform int CarryCoils;
uniform int PrevCoilCount;
uniform vec3 PrevCoilData[ 12 ];
uniform float PrevGuideBz;

layout( location = 0 ) out vec4 outA;
layout( location = 1 ) out vec4 outB;
layout( location = 2 ) out vec4 outC;
layout( location = 3 ) out vec4 outD;

void main()
{
	ivec2 cell = ivec2( gl_FragCoord.xy );
	vec4 A = texelFetch( StateA, cell, 0 );
	vec4 B = texelFetch( StateB, cell, 0 );
	vec4 C = texelFetch( StateC, cell, 0 );
	vec4 D = texelFetch( StateD, cell, 0 );
	Q w    = toPrim( A, B, C, D );

	vec2 x    = cellCentre( cell );
	vec4 clip = clamp( texture( InputTexture, clipUV( cell ) * MaxUV ), 0.0, 1.0 );
	float chi = clamp( w.q3.a, 0.0, 1.0 );
	float luma = dot( clip.rgb, vec3( 0.2126, 0.7152, 0.0722 ) );

	if( CarryCoils == 1 )
	{
		vec2 then = vec2( 0.0 );
		for( int k = 0; k < PrevCoilCount; ++k )
		{
			vec2 d = x - PrevCoilData[ k ].xy;
			then += PrevCoilData[ k ].z * vec2( -d.y, d.x ) / dot( d, d );
		}
		w.q1.yzw += vacuumField( x ) - vec3( then, PrevGuideBz );
	}

	//Feed: inside the ball's footprint the picture relaxes towards the clip.
	w.q3.rgb += FeedFraction * chi * ( clip.rgb - w.q3.rgb );

	//Heating: the clip's light (through Feed) and the audio, in the ball.
	w.q1.x += ( ClipHeatRate * luma + AudioHeatRate ) * chi;

	//Fuel: a steady gas puff. Below the ignition profile, density and
	//pressure are topped back up towards it by FuelFraction this frame; the
	//new gas arrives at rest (momentum unchanged) carrying the clip's colour
	//and marked as ball. It never takes anything away.
	if( FuelFraction > 0.0 )
	{
		float r = length( x - BallCentre );
		float profile = ProfileKind == 0 ? exp( -( r * r ) / ( BallRadius * BallRadius ) )
		                                 : 0.5 - 0.5 * tanh( clamp( ( r - BallRadius ) / ( 0.75 * Dx ), -20.0, 20.0 ) );
		float rhoT = BgDensity + ( 1.0 - BgDensity ) * profile;
		float pT   = BgPressure + ( BallPressure - BgPressure ) * profile;
		float add  = max( rhoT - w.q0.x, 0.0 ) * FuelFraction;
		if( add > 0.0 )
		{
			float rn  = w.q0.x + add;
			w.q0.yzw *= w.q0.x / rn;
			w.q3      = mix( w.q3, vec4( clip.rgb, 1.0 ), add / rn );
			w.q0.x    = rn;
		}
		w.q1.x += max( pT - w.q1.x, 0.0 ) * FuelFraction;
	}

	//A pellet: cold, dense, at rest. Momentum and pressure are unchanged, so
	//the ball has to share its motion with it and it arrives cold.
	if( PelletCount > 0 )
	{
		vec2 d   = x - PelletCentre;
		float add = float( PelletCount ) * PelletDensity * exp( -dot( d, d ) / ( PelletRadius * PelletRadius ) );
		float r   = w.q0.x + add;
		w.q0.yzw *= w.q0.x / r;
		w.q3      = mix( w.q3, vec4( clip.rgb, 1.0 ), add / r );
		w.q0.x    = r;
	}

	toCons( w, A, B, C, D );
	outA = A;
	outB = B;
	outC = vec4( C.x, signalSpeeds( w ).x, C.z, 0.0 );
	outD = D;
}
`;

const EMISSION_MAIN = `
uniform sampler2D StateA;
uniform sampler2D StateB;
uniform sampler2D StateC;
uniform sampler2D StateD;
uniform int View;
uniform float Tint;
uniform int RampKind;       //0 Hot, 1 Aurora
uniform float TempReference;//the temperature at the middle of the ramp
uniform float FieldReference;

out vec4 fragColor;

//Deep red -> violet -> white.
vec3 hotRamp( float s )
{
	vec3 a = vec3( 0.55, 0.03, 0.0 );
	vec3 b = vec3( 0.55, 0.15, 1.0 );
	vec3 c = vec3( 1.0 );
	return s < 0.5 ? mix( a, b, s * 2.0 ) : mix( b, c, s * 2.0 - 1.0 );
}

//The oxygen green line at 557.7 nm, dark when cool, white-green when hot.
vec3 auroraRamp( float s )
{
	vec3 a = vec3( 0.0, 0.22, 0.04 );
	vec3 b = vec3( 0.33, 1.0, 0.18 );
	vec3 c = vec3( 0.86, 1.0, 0.82 );
	return s < 0.5 ? mix( a, b, s * 2.0 ) : mix( b, c, s * 2.0 - 1.0 );
}

vec3 diverging( float s )//-1..1: blue, black, red
{
	return s < 0.0 ? vec3( 0.2, 0.45, 1.0 ) * -s : vec3( 1.0, 0.35, 0.15 ) * s;
}

void main()
{
	//The emission covers the frame, not the margin.
	ivec2 cell = ivec2( gl_FragCoord.xy ) + FrameOrigin;
	Q w = fetchPrim( StateA, StateB, StateC, StateD, cell );
	float r = w.q0.x;
	float p = w.q1.x;
	vec3 b  = w.q1.yzw;
	float T = p / r;

	if( View == 0 )
	{
		//Optically thin bremsstrahlung, integrated along the axis: rho^2 sqrt(T).
		float e  = r * r * sqrt( max( T, 0.0 ) );
		float s  = clamp( 0.5 + 0.25 * log2( max( T, 1e-6 ) / TempReference ), 0.0, 1.0 );
		vec3 ramp = RampKind == 1 ? auroraRamp( s ) : hotRamp( s );
		vec3 col  = mix( max( w.q3.rgb, 0.0 ), ramp, Tint );
		fragColor = vec4( col * e, e );
		return;
	}

	float v = 0.0;
	vec3 col;
	if( View == 1 )
		col = vec3( clamp( ( log( r ) / log( 10.0 ) + 2.0 ) / 2.5, 0.0, 1.0 ) );
	else if( View == 2 )
		col = vec3( clamp( ( log( p ) / log( 10.0 ) + 4.0 ) / 4.5, 0.0, 1.0 ) );
	else if( View == 3 )
		col = vec3( clamp( length( b ) / ( 2.0 * FieldReference ), 0.0, 1.0 ) );
	else if( View == 4 )
		col = diverging( clamp( log( 2.0 * p / max( dot( b, b ), 1e-12 ) ) / log( 10.0 ) / 3.0, -1.0, 1.0 ) );
	else if( View == 5 )
		col = vec3( clamp( length( w.q0.yzw ) / 2.0, 0.0, 1.0 ) );
	else
	{
		//div B dx / |B|, the same measure \`cttest --divb\` asserts, times 10.
		vec3 bxp = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 1, 0 ) ).q1.yzw;
		vec3 bxm = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( -1, 0 ) ).q1.yzw;
		vec3 byp = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 0, 1 ) ).q1.yzw;
		vec3 bym = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 0, -1 ) ).q1.yzw;
		float div = 0.5 * ( ( bxp.x - bxm.x ) + ( byp.y - bym.y ) ) / max( length( b ), 1e-6 );
		col = diverging( clamp( 10.0 * div, -1.0, 1.0 ) );
	}
	fragColor = vec4( col, v );
}
`;

const DOWNSAMPLE_MAIN = `
uniform sampler2D Source;
out vec4 fragColor;
void main()
{
	ivec2 o = ivec2( gl_FragCoord.xy ) * 2;
	fragColor = 0.25 * ( texelFetch( Source, o, 0 ) + texelFetch( Source, o + ivec2( 1, 0 ), 0 )
	                     + texelFetch( Source, o + ivec2( 0, 1 ), 0 ) + texelFetch( Source, o + ivec2( 1, 1 ), 0 ) );
}
`;

const BLUR_MAIN = `
uniform sampler2D Source;
uniform ivec2 Size;
uniform ivec2 Direction;
uniform int Taps;
uniform float Weights[ 129 ];//Weights[ k ] for |offset| = k, normalised over -Taps..Taps
out vec4 fragColor;

int mirrorIndex( int i, int n )
{
	if( i < 0 )
		i = -1 - i;
	if( i >= n )
		i = 2 * n - 1 - i;
	return clamp( i, 0, n - 1 );
}

void main()
{
	ivec2 cell = ivec2( gl_FragCoord.xy );
	vec4 acc   = Weights[ 0 ] * texelFetch( Source, cell, 0 );
	for( int k = 1; k <= Taps; ++k )
	{
		ivec2 a = cell + k * Direction;
		ivec2 b = cell - k * Direction;
		a = ivec2( mirrorIndex( a.x, Size.x ), mirrorIndex( a.y, Size.y ) );
		b = ivec2( mirrorIndex( b.x, Size.x ), mirrorIndex( b.y, Size.y ) );
		acc += Weights[ k ] * ( texelFetch( Source, a, 0 ) + texelFetch( Source, b, 0 ) );
	}
	fragColor = acc;
}
`;

const POISSON_COMMON = `
uniform sampler2D Solution;//A on this level
uniform sampler2D Rhs;     //f on this level
uniform ivec2 Size;
uniform float Spacing;     //this level's cell size
uniform int VacuumGhosts;  //1 on the finest level: A off the grid is the coils'
uniform ivec2 GridSize;
uniform float Dx;
uniform int CoilCount;
uniform vec3 CoilData[ 12 ];

float vacuumPotential( vec2 x )
{
	float a = 0.0;
	for( int k = 0; k < CoilCount; ++k )
	{
		vec2 d = x - CoilData[ k ].xy;
		a -= CoilData[ k ].z * 0.5 * log( dot( d, d ) );
	}
	return a;
}

float solutionAt( ivec2 c )
{
	if( c.x >= 0 && c.y >= 0 && c.x < Size.x && c.y < Size.y )
		return texelFetch( Solution, c, 0 ).r;
	if( VacuumGhosts == 1 )
		return vacuumPotential( ( vec2( c ) + 0.5 ) * Spacing );
	//A correction is zero on the boundary; the ghost is its mirror.
	ivec2 m = clamp( c, ivec2( 0 ), Size - 1 );
	return -texelFetch( Solution, m, 0 ).r;
}

float neighbours( ivec2 c )
{
	return solutionAt( c + ivec2( 1, 0 ) ) + solutionAt( c + ivec2( -1, 0 ) ) + solutionAt( c + ivec2( 0, 1 ) )
	       + solutionAt( c + ivec2( 0, -1 ) );
}
`;

const RESIDUAL_MAIN = `
uniform int Mode;//0: the right-hand side, -J_z, from the state's field; 1: f - del^2 A
uniform sampler2D StateB;
out vec4 fragColor;

vec3 fieldAt( ivec2 c )
{
	return texelFetch( StateB, clamp( c, ivec2( 0 ), Size - 1 ), 0 ).yzw;
}

void main()
{
	ivec2 c = ivec2( gl_FragCoord.xy );
	if( Mode == 0 )
	{
		float h  = 0.5 / Spacing;
		float jz = ( fieldAt( c + ivec2( 1, 0 ) ).y - fieldAt( c + ivec2( -1, 0 ) ).y ) * h
		           - ( fieldAt( c + ivec2( 0, 1 ) ).x - fieldAt( c + ivec2( 0, -1 ) ).x ) * h;
		fragColor = vec4( -jz, 0.0, 0.0, 0.0 );
		return;
	}
	float a   = texelFetch( Solution, c, 0 ).r;
	float lap = ( neighbours( c ) - 4.0 * a ) / ( Spacing * Spacing );
	fragColor = vec4( texelFetch( Rhs, c, 0 ).r - lap, 0.0, 0.0, 0.0 );
}
`;

const SMOOTH_MAIN = `
out vec4 fragColor;
void main()
{
	ivec2 c   = ivec2( gl_FragCoord.xy );
	float a   = texelFetch( Solution, c, 0 ).r;
	float jac = 0.25 * ( neighbours( c ) - Spacing * Spacing * texelFetch( Rhs, c, 0 ).r );
	fragColor = vec4( mix( a, jac, 0.8 ), 0.0, 0.0, 0.0 );
}
`;

const PROLONG_MAIN = `
uniform sampler2D Coarse;
out vec4 fragColor;
void main()
{
	ivec2 c = ivec2( gl_FragCoord.xy );
	vec2 uv = ( vec2( c ) + 0.5 ) / vec2( Size );
	fragColor = vec4( texelFetch( Solution, c, 0 ).r + texture( Coarse, uv ).r, 0.0, 0.0, 0.0 );
}
`;

const COMPOSITE_MAIN = `
uniform sampler2D InputTexture;
uniform sampler2D Emission;
uniform sampler2D Glow0;
uniform sampler2D Glow1;
uniform sampler2D Glow2;
uniform sampler2D Potential;
uniform vec2 MaxUV;
uniform ivec2 OutputSize;
uniform int ExactInput;    //1: the host texture is exactly the output's size
uniform int View;
uniform float Gain;        //2^EV / the reference emission
uniform float GlowAmount;  //the fraction moved into the glare
uniform float CoreWeight;  //1 - GlowAmount: what the glare took from the core
uniform vec3 GlowShare;
uniform float LineAmount;
uniform float LineSpacing; //A_z between two drawn field lines; 0 = none
uniform vec4 PotentialMap; //frame uv -> the grid's: uv * xy + zw (the potential covers the margin)
uniform float MixAmount;

in vec2 uv;
out vec4 fragColor;

//The fleet's shoulder: linear to the knee, then an exponential approach to 1.
//Not a clamp, whose corner would flatten the core of the orb into a disc.
float shoulder( float v )
{
	const float knee = 0.8;
	if( v <= knee )
		return v;
	float k = 1.0 - knee;
	return knee + k * ( 1.0 - exp( -( v - knee ) / k ) );
}

void main()
{
	vec4 source;
	if( ExactInput == 1 )
		source = texelFetch( InputTexture, ivec2( gl_FragCoord.xy ), 0 );
	else
		source = texture( InputTexture, uv * MaxUV );

	//Mix 0 is the identity, bit for bit.
	if( MixAmount <= 0.0 )
	{
		fragColor = source;
		return;
	}

	vec4 e = texture( Emission, uv );
	vec3 col;
	if( View == 0 )
	{
		vec3 light = CoreWeight * e.rgb;
		light += GlowAmount * ( GlowShare.x * texture( Glow0, uv ).rgb + GlowShare.y * texture( Glow1, uv ).rgb
		                        + GlowShare.z * texture( Glow2, uv ).rgb );

		//Field lines: contours of A_z, one pixel wide wherever they are, lit
		//by the plasma sitting on them.
		if( LineAmount > 0.0 && LineSpacing > 0.0 )
		{
			float f  = texture( Potential, uv * PotentialMap.xy + PotentialMap.zw ).r / LineSpacing;
			float fw = max( fwidth( f ), 1e-6 );
			float d  = abs( fract( f + 0.5 ) - 0.5 ) / fw;
			//Where the lines are closer than about two pixels they cannot be
			//drawn as lines -- every pixel is within a pixel of one and the
			//whole region just lights up (it did, beside the coils, and it
			//made Field Line Count do nothing). They fade out there instead.
			float resolved = clamp( 1.5 - 3.0 * fw, 0.0, 1.0 );
			light += LineAmount * 3.0 * resolved * exp( -d * d ) * e.rgb;
		}

		light *= Gain;
		col = vec3( shoulder( light.r ), shoulder( light.g ), shoulder( light.b ) );
	}
	else
		col = e.rgb;

	fragColor = mix( source, vec4( col, 1.0 ), MixAmount );
}
`;


/// `SourceFor()`: every program as the plugin assembles it. The vertex shader
/// is `kVersion + kVertex` for all of them. check_shaders.py compares each
/// row with the matching `case` in Shaders.cpp.
const VERTEX_SOURCE = VERSION + VERTEX;
const PROGRAM_PIECES = {
  reduce: [VERSION, REDUCE_MAIN],
  clock: [VERSION, CLOCK_MAIN],
  predict: [VERSION, COMMON, RIEMANN, FORCES, PREDICT_MAIN],
  flux: [VERSION, COMMON, RIEMANN, FLUX_MAIN],
  update: [VERSION, COMMON, RIEMANN, FORCES, UPDATE_MAIN],
  ignite: [VERSION, COMMON, IGNITE_MAIN],
  sources: [VERSION, COMMON, SOURCES_MAIN],
  emission: [VERSION, COMMON, EMISSION_MAIN],
  downsample: [VERSION, DOWNSAMPLE_MAIN],
  blur: [VERSION, BLUR_MAIN],
  residual: [VERSION, POISSON_COMMON, RESIDUAL_MAIN],
  smooth: [VERSION, POISSON_COMMON, SMOOTH_MAIN],
  prolong: [VERSION, POISSON_COMMON, PROLONG_MAIN],
  composite: [VERSION, COMPOSITE_MAIN],
};

//===========================================================================
// Controls.h / Controls.cpp, ported.
//===========================================================================

const f32 = Math.fround;
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
const clamp01 = (v) => clamp(v, 0, 1);
const geometric = (v, lo, hi) => f32(lo * Math.pow(hi / lo, clamp01(v)));
const linear = (v, lo, hi) => f32(lo + (hi - lo) * clamp01(v));

const K_SPEED_LOW = 0.02;
const K_SPEED_HIGH = 2.0;

const ballSizeFromParam = (v) => linear(v, 0.04, 0.4);
const temperatureFromParam = (v) => geometric(v, 0.1, 50.0);
const feedFromParam = (v) => clamp01(v);
const clipHeatsFromParam = (v) => clamp01(v);
const fuelFromParam = (v) => linear(v, 0.0, 2.0);
const fieldFromParam = (v) => geometric(v, 0.25, 4.0);
const guideFieldFromParam = (v) => linear(v, 0.0, 2.0);
const coilRadiusFromParam = (v) => linear(v, 1.1, 3.0);
const coilSpinFromParam = (v) => linear(v, -1.0, 1.0);
const curvatureFromParam = (v) => { const c = clamp01(v); return f32(3.0 * c * c); };
const speedFromParam = (v) => (v <= 0 ? 0 : geometric(v, K_SPEED_LOW, K_SPEED_HIGH));
function paramFromSpeed(speed) {
  if (speed <= 0) return 0;
  const clamped = clamp(speed, K_SPEED_LOW, K_SPEED_HIGH);
  return f32(Math.log(clamped / K_SPEED_LOW) / Math.log(K_SPEED_HIGH / K_SPEED_LOW));
}
const resistivityFromParam = (v) => (v <= 0 ? 0 : geometric(v, 1.0e-5, 1.0e-2));
const coolingFromParam = (v) => { const c = clamp01(v); return f32(8.0 * c * c * c); };
const driveFromParam = (v) => { const c = clamp01(v); return f32(4.0 * c * c); };
const driveScaleFromParam = (v) => geometric(v, 0.05, 0.6);
const exposureFromParam = (v) => linear(v, -4.0, 8.0);
const glowFromParam = (v) => linear(v, 0.0, 0.9);
const lineCountFromParam = (v) => clamp(Math.round(linear(v, 4.0, 48.0)), 4, 48);

/// Controls.h: the model's constants. None of these is a control.
const POLE_COUNTS = [0, 4, 6, 8, 12];
const DETAIL_CELLS = [128, 256, 512, 1024];
const K_GAMMA = f32(5.0 / 3.0);
const K_BACKGROUND_DENSITY = f32(0.1);
const K_BACKGROUND_PRESSURE = f32(0.005);
const K_DENSITY_FLOOR = f32(1.0e-4);
const K_PRESSURE_FLOOR = f32(1.0e-6);
const K_COURANT = f32(0.8);
const K_MAX_SUBSTEPS = 48;
const K_MARGIN_FRACTION = f32(0.1);
const K_SPONGE_EFOLDS = f32(6.0);
const K_SPONGE_POWER = f32(2.0);
const K_ENTROPY_SWITCH = f32(0.02);
const K_GLM_ALPHA = f32(0.4);
const K_PELLET_DENSITY = f32(2.0);
const K_PELLET_RADIUS = f32(0.035);
const K_QUENCH_TIME = f32(0.15);
const K_DRIVE_TIME = f32(0.6);
const K_DRIVE_MODES = 12;
const K_DRIVE_WINDOW = f32(1.5);
const K_EMISSION_REFERENCE = f32(0.5);
const BOUNDARY_OPEN = 0;
const BOUNDARY_WALL = 1;

/// Shaders.h.
const K_GLOW_STAGES = 3;
const K_GLOW_SIGMA = [f32(0.012), f32(0.045), f32(0.15)];
const K_GLOW_SHARE = [f32(0.5), f32(0.33), f32(0.17)];
const K_GLOW_CELLS = 256;
const K_MAX_BLUR_TAPS = 128;

/// Containment.cpp's constants.
const K_REDUCE_BLOCK = 4;
const K_MAX_FRAME_DELTA = 0.1;
const K_MAX_COILS = 12;
const K_PI = 3.14159265358979323846;

/// The one departure from the constructor's defaults: Detail 128 rather than
/// 256. See the note at the top and the page's disclosure.
const DEMO_DETAIL = 0;

//===========================================================================
// Physics.cpp, ported: the CPU half.
//===========================================================================

/// `VacuumField` (//= mirrored in the shader's vacuumField()).
function vacuumField(coils, x, y) {
  let bx = 0;
  let by = 0;
  for (let k = 0; k < coils.count; k += 1) {
    const dx = x - coils.x[k];
    const dy = y - coils.y[k];
    const d2 = dx * dx + dy * dy;
    bx += -coils.current[k] * dy / d2;
    by += coils.current[k] * dx / d2;
  }
  return [bx, by, coils.guide];
}

/// `VacuumPotential` (//= mirrored in the shader's vacuumPotential()).
function vacuumPotential(coils, x, y) {
  let a = 0;
  for (let k = 0; k < coils.count; k += 1) {
    const dx = x - coils.x[k];
    const dy = y - coils.y[k];
    a -= coils.current[k] * 0.5 * Math.log(dx * dx + dy * dy);
  }
  return a;
}

function emptyCoils() {
  return { count: 0, x: new Float64Array(K_MAX_COILS), y: new Float64Array(K_MAX_COILS),
    current: new Float64Array(K_MAX_COILS), guide: 0 };
}

/// `MakeCoils`: alternating line currents on a circle, normalised so |B| at
/// radius 0.5 in the first gap is `field` (before `strength`).
function makeCoils(poles, radius, angle, field, guide, strength, cx, cy) {
  const coils = emptyCoils();
  coils.count = clamp(poles, 0, K_MAX_COILS);
  coils.guide = guide * strength;
  if (coils.count === 0) return coils;

  for (let k = 0; k < coils.count; k += 1) {
    const phi = angle + 2.0 * K_PI * k / coils.count;
    coils.x[k] = cx + radius * Math.cos(phi);
    coils.y[k] = cy + radius * Math.sin(phi);
    coils.current[k] = k % 2 === 0 ? 1.0 : -1.0;
  }
  const gap = angle + K_PI / coils.count;
  const [bx, by] = vacuumField(coils, cx + 0.5 * Math.cos(gap), cy + 0.5 * Math.sin(gap));
  const here = Math.sqrt(bx * bx + by * by);
  const scale = here > 0 ? field * strength / here : 0;
  for (let k = 0; k < coils.count; k += 1) coils.current[k] *= scale;
  return coils;
}

function copyCoils(c) {
  return { count: c.count, x: Float64Array.from(c.x), y: Float64Array.from(c.y),
    current: Float64Array.from(c.current), guide: c.guide };
}

/// `ChooseGrid`: square cells, `cells` on the frame's short side, the long side
/// a multiple of 8, and Open's margin all round.
function chooseGrid(width, height, cells, margin = 0) {
  const landscape = width >= height;
  const aspect = landscape ? width / Math.max(height, 1) : height / Math.max(width, 1);
  const shortSide = cells;
  const longSide = Math.max(8, 8 * Math.round(shortSide * aspect / 8.0));
  const g = {};
  g.fx = landscape ? longSide : shortSide;
  g.fy = landscape ? shortSide : longSide;
  g.ox = g.oy = Math.max(margin, 0);
  g.nx = g.fx + 2 * g.ox;
  g.ny = g.fy + 2 * g.oy;
  g.dx = 1.0 / g.fy;
  g.lx = g.nx * g.dx;
  g.ly = g.ny * g.dx;
  return g;
}

/// `Pcg`: PCG-XSH-RR on a 64-bit state, in BigInt so the stirring wanders the
/// way the plugin's does from the same seed.
const MASK64 = (1n << 64n) - 1n;
class Pcg {
  constructor(seed = 0x853c49e6748fea9bn) {
    this.increment = 0xda3e39cb94b95bdbn;
    this.state = 0n;
    this.next();
    this.state = (this.state + seed) & MASK64;
    this.next();
  }

  next() {
    const old = this.state;
    this.state = (old * 6364136223846793005n + this.increment) & MASK64;
    const xorshifted = Number((((old >> 18n) ^ old) >> 27n) & 0xffffffffn);
    const rot = Number(old >> 59n);
    return ((xorshifted >>> rot) | (xorshifted << ((-rot) & 31))) >>> 0;
  }

  uniform() {
    return this.next() / 4294967296.0;
  }

  normal() {
    const u1 = Math.max(this.uniform(), 1e-12);
    const u2 = this.uniform();
    return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * K_PI * u2);
  }
}

/// `Drive`: twelve Fourier modes of a stream function, each amplitude an
/// Ornstein-Uhlenbeck process with correlation time kDriveTime.
class Drive {
  constructor() {
    this.random = new Pcg();
    this.modes = [];
    this.state = [];
    this.angle = [];
    this.band = [];
  }

  reset() {
    this.random = new Pcg(0x5eedn);
    this.modes = [];
    this.state = [];
    this.angle = [];
    this.band = [];
  }

  advance(dt, scale, strength) {
    if (this.state.length === 0) {
      for (let i = 0; i < 2 * K_DRIVE_MODES; i += 1) this.state.push(this.random.normal());
      for (let m = 0; m < K_DRIVE_MODES; m += 1) {
        this.angle.push(2.0 * K_PI * this.random.uniform());
        this.band.push(0.75 + 0.5 * this.random.uniform());
      }
    }
    const decay = Math.exp(-Math.max(dt, 0.0) / K_DRIVE_TIME);
    const kick = Math.sqrt(Math.max(0.0, 1.0 - decay * decay));
    for (let i = 0; i < this.state.length; i += 1) this.state[i] = this.state[i] * decay + kick * this.random.normal();

    this.modes = [];
    for (let m = 0; m < K_DRIVE_MODES; m += 1) {
      const k = 2.0 * K_PI / Math.max(scale, 1e-3) * this.band[m];
      const amp = strength / (k * Math.sqrt(K_DRIVE_MODES));
      this.modes.push({
        kx: f32(k * Math.cos(this.angle[m])),
        ky: f32(k * Math.sin(this.angle[m])),
        a: f32(amp * this.state[2 * m]),
        b: f32(amp * this.state[2 * m + 1]),
      });
    }
  }
}

//===========================================================================
// Presets.h, ported: the rows, and which parameter each column drives
// (`kPresetTarget`). Standard columns hold the host's 0..1; option and
// boolean columns hold their element value.
//===========================================================================

const PRESET_COLUMNS = [
  'ballSize', 'temperature', 'profile', 'feed', 'clipHeats', 'fuel', 'field', 'guideField',
  'poles', 'coilRadius', 'coilSpin', 'curvature', 'quench', 'boundary', 'resistivity', 'cooling',
  'drive', 'driveScale', 'exposure', 'tint', 'ramp', 'glow', 'fieldLines', 'lineCount',
];

// clang-format off
//                                size   temp   prof feed  clip  fuel  field guide poles radius spin  curv   q  bnd res cool drive dscale expo    tint ramp glow   lines count
const PRESETS = [
  { name: 'Clip Orb',           v: [0.389, 0.315, 0, 0.1, 0.0, 0.15, 0.25, 0.7, 2, 0.105, 0.55, 0.316, 0, 1, 0, 0, 0.274, 0.442, 0.4167, 0.0, 0, 0.333, 0.0, 0.273] },
  { name: 'Green Orb',          v: [0.389, 0.315, 0, 0.1, 0.0, 0.15, 0.25, 0.7, 2, 0.105, 0.62, 0.36, 0, 1, 0, 0, 0.33, 0.4, 0.45, 0.9, 1, 0.62, 0.35, 0.25] },
  { name: 'Guide Field Bubble', v: [0.306, 0.4, 1, 0.1, 0.0, 0.0, 0.5, 0.5, 0, 0.105, 0.5, 0.0, 0, 1, 0, 0, 0.0, 0.442, 0.4167, 0.5, 0, 0.4, 0.0, 0.273] },
  { name: 'Cusp Leak',          v: [0.222, 0.482, 0, 0.1, 0.0, 0.25, 0.5, 0.0, 2, 0.105, 0.5, 0.0, 0, 0, 0, 0, 0.0, 0.442, 0.4167, 0.4, 0, 0.45, 0.3, 0.273] },
  { name: 'Rayleigh-Taylor',    v: [0.389, 0.315, 1, 0.1, 0.0, 0.15, 0.25, 0.7, 0, 0.105, 0.5, 0.7, 0, 1, 0, 0, 0.1, 0.442, 0.4167, 0.3, 0, 0.4, 0.0, 0.273] },
  { name: 'Quench Fireball',    v: [0.306, 0.741, 1, 0.0, 0.0, 0.0, 0.5, 0.7, 2, 0.105, 0.5, 0.0, 1, 0, 0, 0, 0.0, 0.442, 0.4167, 1.0, 0, 0.6, 0.0, 0.273] },
];
// clang-format on

//===========================================================================
// Buffers.
//===========================================================================

/// `StateBuffer`: four RGBA32F targets on one framebuffer, nearest, clamped,
/// cleared to zero when made.
class StateBuffer {
  constructor(gl) {
    this.gl = gl;
    this.width = 0;
    this.height = 0;
    this.textures = [];
    this.fbo = null;
  }

  ensure(w, h) {
    const gl = this.gl;
    if (this.fbo && w === this.width && h === this.height) return;
    this.dispose();
    for (let i = 0; i < 4; i += 1) {
      const t = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA32F, w, h);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      this.textures.push(t);
    }
    gl.bindTexture(gl.TEXTURE_2D, null);
    this.fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    const attachments = [gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2, gl.COLOR_ATTACHMENT3];
    attachments.forEach((a, i) => gl.framebufferTexture2D(gl.FRAMEBUFFER, a, gl.TEXTURE_2D, this.textures[i], 0));
    gl.drawBuffers(attachments);
    const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (status === gl.FRAMEBUFFER_COMPLETE) {
      gl.viewport(0, 0, w, h);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    if (status !== gl.FRAMEBUFFER_COMPLETE) {
      this.dispose();
      throw new GLError(`the plasma's state buffer is incomplete (0x${status.toString(16)}) at ${w}x${h}`);
    }
    this.width = w;
    this.height = h;
  }

  dispose() {
    const gl = this.gl;
    if (this.fbo) gl.deleteFramebuffer(this.fbo);
    for (const t of this.textures) gl.deleteTexture(t);
    this.fbo = null;
    this.textures = [];
    this.width = 0;
    this.height = 0;
  }
}

//===========================================================================
// What the stats line under the canvas reports. Filled by the renderer.
//===========================================================================

const telemetry = {
  grid: '', substeps: 0, target: 0, done: 0, capped: 0, frames: 0, floors: 0, entropy: 0, cells: 0,
};

//===========================================================================
// The plugin's frame, in the order ProcessOpenGL runs it.
//===========================================================================

class ContainmentRenderer {
  constructor(gl, quad) {
    this.gl = gl;
    this.quad = quad;

    // The emission, the glow and the potential are RGBA32F read with a LINEAR
    // filter, as in the plugin. Without OES_texture_float_linear they would be
    // incomplete and sample as black: an orb that is not there, rather than an
    // error. So its absence stops the page.
    if (!gl.getExtension('OES_texture_float_linear')) {
      throw new GLError('OES_texture_float_linear is missing. The plasma\'s light, its glow and the field lines\' potential are 32-bit float textures read with a linear filter, as in the plugin; without the extension they would sample as black.');
    }
    // The update pass reads fifteen textures at once and every state pass
    // writes four. Both are within WebGL2's guaranteed minimums (16 and 4),
    // but say so plainly if a browser offers less.
    if (gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS) < 15) {
      throw new GLError('This GPU offers fewer than 15 texture units to a fragment shader. The plugin\'s update pass reads fifteen.');
    }
    if (gl.getParameter(gl.MAX_DRAW_BUFFERS) < 4) {
      throw new GLError('This GPU cannot write four render targets at once. The plugin\'s state is four RGBA32F textures written together.');
    }

    this.programs = {};
    for (const [name, pieces] of Object.entries(PROGRAM_PIECES)) {
      this.programs[name] = new Program(gl, VERTEX_SOURCE, pieces.join(''), name);
    }

    this.state = [new StateBuffer(gl), new StateBuffer(gl)];
    this.star = new StateBuffer(gl);
    this.fluxX = new StateBuffer(gl);
    this.fluxY = new StateBuffer(gl);
    this.current = 0;
    this.reduction = [];
    this.clock = [new PassBuffer(gl, { filter: 'nearest' }), new PassBuffer(gl, { filter: 'nearest' })];
    this.clockIndex = 0;

    this.emission = new PassBuffer(gl, { filter: 'linear' });
    this.glowReduce = [];
    this.glowTemp = new PassBuffer(gl, { filter: 'nearest' });
    this.glow = [0, 1, 2].map(() => new PassBuffer(gl, { filter: 'linear' }));
    this.glowWidth = 0;
    this.glowHeight = 0;
    this.glowFraction = 0;

    this.levels = [];
    this.havePotential = false;

    this.grid = { nx: 0, ny: 0, fx: 0, fy: 0, ox: 0, oy: 0, dx: 0, lx: 0, ly: 0 };
    this.coils = emptyCoils();
    this.previousCoils = emptyCoils();
    this.carryCoils = false;
    this.drive = new Drive();
    this.spinAngle = 0;
    this.coilStrength = 1;
    this.frameSubsteps = 0;
    this.lastSpeed = 0;
    this.speedKnown = false;
    this.clockPending = false;
    this.lastTarget = 0;
    this.lastEntropyCells = 0;

    this.lastNow = -1;
    this.simTime = 0;

    this.ignitePending = true;
    this.pelletPresses = 0;
    this.lastPreset = null;
    this.params = null;

    this.readback = new Float32Array(8);
  }

  //-------------------------------------------------------------------------
  // `P()`: the operator's value, or the preset row's where the Preset
  // dropdown is on anything but Custom.
  //-------------------------------------------------------------------------
  P(id) {
    const preset = clamp(Math.round(this.params.get('preset')), 0, PRESETS.length);
    if (preset > 0 && id !== 'preset') {
      const column = PRESET_COLUMNS.indexOf(id);
      if (column >= 0) return PRESETS[preset - 1].v[column];
    }
    return this.params.get(id);
  }

  option(id, count) {
    return clamp(Math.round(this.P(id)), 0, count - 1);
  }

  //-------------------------------------------------------------------------
  // The events. The plugin counts rising edges in SetFloatParameter; here the
  // button is a boolean, taken and released on the frame it is seen. A new
  // Preset re-ignites, as SetFloatParameter does.
  //-------------------------------------------------------------------------
  takeEvents(params) {
    if (params.get('ignite') > 0.5) { this.ignitePending = true; params.set('ignite', 0); }
    if (params.get('pellet') > 0.5) { this.pelletPresses += 1; params.set('pellet', 0); }
    const preset = Math.round(params.get('preset'));
    if (this.lastPreset !== null && preset !== this.lastPreset) this.ignitePending = true;
    this.lastPreset = preset;
  }

  //-------------------------------------------------------------------------
  // Uniform helpers: the kit's Program.set() is float-only; the plugin's
  // integer uniforms go through these, as its uniform1i / uniform2i do.
  //-------------------------------------------------------------------------
  i1(p, name, v) { p.setInt(name, v); }

  i2(p, name, a, b) {
    const loc = p.location(name);
    if (loc !== null) this.gl.uniform2i(loc, a, b);
  }

  f1(p, name, v) { p.set(name, v); }

  f2(p, name, a, b) { p.set(name, a, b); }

  samplers(p, names) { names.forEach((n, unit) => p.setSampler(n, unit)); }

  /// `BoundTextures`: bind to units 0..n-1, draw, and unbind them all again.
  withTextures(textures, draw) {
    const gl = this.gl;
    textures.forEach((t, unit) => bindTexture(gl, unit, t));
    gl.activeTexture(gl.TEXTURE0);
    draw();
    for (let unit = textures.length - 1; unit >= 0; unit -= 1) bindTexture(gl, unit, null);
    gl.activeTexture(gl.TEXTURE0);
  }

  target(buffer, width, height) {
    const gl = this.gl;
    gl.bindFramebuffer(gl.FRAMEBUFFER, buffer.fbo);
    gl.viewport(0, 0, width ?? buffer.width, height ?? buffer.height);
  }

  //-------------------------------------------------------------------------
  // `CurrentModel()`.
  //-------------------------------------------------------------------------
  currentModel() {
    const g = this.grid;
    const m = {};
    m.field = fieldFromParam(this.P('field'));
    m.guide = guideFieldFromParam(this.P('guideField')) * m.field;
    m.poles = POLE_COUNTS[this.option('poles', POLE_COUNTS.length)];
    m.spin = coilSpinFromParam(this.P('coilSpin'));
    m.curvature = curvatureFromParam(this.P('curvature'));
    m.boundary = this.option('boundary', 2);
    m.eta = resistivityFromParam(this.P('resistivity'));
    m.cooling = coolingFromParam(this.P('cooling'));
    m.drive = driveFromParam(this.P('drive'));
    m.driveScale = driveScaleFromParam(this.P('driveScale'));
    const halfDiagonal = 0.5 * Math.sqrt(g.lx * g.lx + g.ly * g.ly);
    m.coilRadius = Math.max(coilRadiusFromParam(this.P('coilRadius')), 1.05 * halfDiagonal);
    m.ballX = (g.ox + clamp01(this.P('ballX')) * g.fx) * g.dx;
    m.ballY = (g.oy + clamp01(this.P('ballY')) * g.fy) * g.dx;
    m.ballRadius = ballSizeFromParam(this.P('ballSize'));
    m.pressure = 0.5 * temperatureFromParam(this.P('temperature'));
    return m;
  }

  //-------------------------------------------------------------------------
  // `EnsureBuffers()`.
  //-------------------------------------------------------------------------
  ensureBuffers(wanted) {
    const gl = this.gl;
    const reset = wanted.nx !== this.grid.nx || wanted.ny !== this.grid.ny;
    if (reset) this.ignitePending = true;

    this.state[0].ensure(wanted.nx, wanted.ny);
    this.state[1].ensure(wanted.nx, wanted.ny);
    this.star.ensure(wanted.nx, wanted.ny);
    this.fluxX.ensure(wanted.nx + 1, wanted.ny);
    this.fluxY.ensure(wanted.nx, wanted.ny + 1);

    if (reset || this.reduction.length === 0) {
      for (const b of this.reduction) b.dispose();
      this.reduction = [];
      let w = wanted.nx;
      let h = wanted.ny;
      do {
        w = Math.ceil(w / K_REDUCE_BLOCK);
        h = Math.ceil(h / K_REDUCE_BLOCK);
        this.reduction.push(new PassBuffer(gl, { filter: 'nearest' }).ensure(w, h, gl.RGBA32F));
      } while (w * h > 16);
    }
    for (const b of this.clock) b.ensure(2, 1, gl.RGBA32F);

    this.emission.ensure(wanted.fx, wanted.fy, gl.RGBA32F);
    let gw = wanted.fx;
    let gh = wanted.fy;
    let reductions = 0;
    while (Math.min(gw, gh) > K_GLOW_CELLS && gw % 2 === 0 && gh % 2 === 0) {
      gw /= 2;
      gh /= 2;
      reductions += 1;
    }
    if (this.glowReduce.length !== reductions) {
      for (const b of this.glowReduce) b.dispose();
      this.glowReduce = [];
      for (let i = 0; i < reductions; i += 1) this.glowReduce.push(new PassBuffer(gl, { filter: 'nearest' }));
    }
    {
      let w = wanted.fx;
      let h = wanted.fy;
      for (const b of this.glowReduce) {
        w /= 2;
        h /= 2;
        b.ensure(w, h, gl.RGBA32F);
      }
    }
    this.glowTemp.ensure(gw, gh, gl.RGBA32F);
    for (const b of this.glow) b.ensure(gw, gh, gl.RGBA32F);
    this.glowWidth = gw;
    this.glowHeight = gh;

    if (reset || this.levels.length === 0) {
      for (const level of this.levels) {
        level.solution[0].dispose();
        level.solution[1].dispose();
        level.rhs.dispose();
        level.residual.dispose();
      }
      this.levels = [];
      let w = wanted.nx;
      let h = wanted.ny;
      let spacing = wanted.dx;
      for (;;) {
        this.levels.push({
          solution: [new PassBuffer(gl, { filter: 'linear' }), new PassBuffer(gl, { filter: 'linear' })],
          rhs: new PassBuffer(gl, { filter: 'nearest' }),
          residual: new PassBuffer(gl, { filter: 'nearest' }),
          current: 0, nx: w, ny: h, spacing,
        });
        if (w % 2 !== 0 || h % 2 !== 0 || Math.min(w, h) < 16) break;
        w /= 2;
        h /= 2;
        spacing *= 2.0;
      }
      this.havePotential = false;
    }
    for (const level of this.levels) {
      level.solution[0].ensure(level.nx, level.ny, gl.RGBA32F);
      level.solution[1].ensure(level.nx, level.ny, gl.RGBA32F);
      level.rhs.ensure(level.nx, level.ny, gl.RGBA32F);
      level.residual.ensure(level.nx, level.ny, gl.RGBA32F);
    }

    this.grid = { ...wanted };
  }

  //-------------------------------------------------------------------------
  // `SetStateUniforms()`.
  //-------------------------------------------------------------------------
  setStateUniforms(p) {
    const g = this.grid;
    this.i2(p, 'GridSize', g.nx, g.ny);
    this.f1(p, 'Dx', g.dx);
    this.i2(p, 'FrameOrigin', g.ox, g.oy);
    this.i2(p, 'FrameCells', g.fx, g.fy);
    this.f1(p, 'Gamma', K_GAMMA);
    const m = this.currentModel();
    this.i1(p, 'BoundaryMode', m.boundary);
    this.i1(p, 'UseFloors', 1);
    this.f1(p, 'RhoFloor', K_DENSITY_FLOOR);
    this.f1(p, 'PFloor', K_PRESSURE_FLOOR);
    this.f1(p, 'AmbientDensity', K_BACKGROUND_DENSITY);
    this.f1(p, 'AmbientPressure', K_BACKGROUND_PRESSURE);
    this.i1(p, 'Solver', 0);
    this.i1(p, 'UseGLM', 1);

    p.setArray('CoilData', coilData(this.coils), 3);
    this.i1(p, 'CoilCount', this.coils.count);
    this.f1(p, 'GuideBz', this.coils.guide);

    this.f1(p, 'Curvature', m.curvature);
    this.f1(p, 'CurvatureTemp', m.pressure);
    this.f2(p, 'GravityCentre', m.ballX, m.ballY);
    this.f1(p, 'GravityCore', 0.5 * m.ballRadius);

    const modes = new Float32Array(K_DRIVE_MODES * 4);
    let count = 0;
    if (m.drive > 0) {
      for (const mode of this.drive.modes) {
        modes[count * 4 + 0] = mode.kx;
        modes[count * 4 + 1] = mode.ky;
        modes[count * 4 + 2] = mode.a;
        modes[count * 4 + 3] = mode.b;
        count += 1;
      }
    }
    this.i1(p, 'DriveCount', count);
    this.f1(p, 'DriveWindow', K_DRIVE_WINDOW * m.ballRadius);
    p.setArray('DriveModes', modes, 4);
  }

  drawInto(target) {
    this.target(target);
    this.quad.draw();
  }

  //-------------------------------------------------------------------------
  // `Ignite()`.
  //-------------------------------------------------------------------------
  ignite(input) {
    const m = this.currentModel();
    const p = this.programs.ignite.use();
    this.setStateUniforms(p);
    this.samplers(p, ['InputTexture']);
    this.f2(p, 'MaxUV', 1, 1);
    this.f2(p, 'BallCentre', m.ballX, m.ballY);
    this.f1(p, 'BallRadius', m.ballRadius);
    this.i1(p, 'ProfileKind', this.option('profile', 2));
    this.f1(p, 'BallPressure', m.pressure);
    this.f1(p, 'ClipHeats', clipHeatsFromParam(this.P('clipHeats')));
    this.f1(p, 'BgDensity', K_BACKGROUND_DENSITY);
    this.f1(p, 'BgPressure', K_BACKGROUND_PRESSURE);
    this.withTextures([input], () => this.drawInto(this.state[this.current]));

    this.drive.reset();
    this.speedKnown = false;
    this.havePotential = false;
  }

  //-------------------------------------------------------------------------
  // `Sources()`, with the analyser's silence written in: Level() is 0, so
  // Audio Heat adds nothing.
  //-------------------------------------------------------------------------
  sources(input, hostDt, simDt, pellets) {
    const m = this.currentModel();
    const feed = feedFromParam(this.P('feed'));
    const fraction = f32(1.0 - Math.pow(1.0 - Math.min(feed, f32(0.999)), hostDt * 60.0));
    const clipHeat = f32(2.0 * clipHeatsFromParam(this.P('clipHeats')) * feed * m.pressure * simDt);
    const audioHeat = 0;
    const fuel = f32(1.0 - Math.exp(-fuelFromParam(this.P('fuel')) * simDt));
    if (fraction <= 0 && clipHeat <= 0 && audioHeat <= 0 && pellets === 0 && fuel <= 0 && !this.carryCoils) return;

    const p = this.programs.sources.use();
    this.setStateUniforms(p);
    this.samplers(p, ['StateA', 'StateB', 'StateC', 'StateD', 'InputTexture']);
    this.f2(p, 'MaxUV', 1, 1);
    this.f1(p, 'FeedFraction', fraction);
    this.f1(p, 'ClipHeatRate', clipHeat);
    this.f1(p, 'AudioHeatRate', audioHeat);
    this.i1(p, 'PelletCount', pellets);
    this.f1(p, 'FuelFraction', fuel);
    this.f2(p, 'BallCentre', m.ballX, m.ballY);
    this.f1(p, 'BallRadius', m.ballRadius);
    this.i1(p, 'ProfileKind', this.option('profile', 2));
    this.f1(p, 'BallPressure', m.pressure);
    this.f1(p, 'BgDensity', K_BACKGROUND_DENSITY);
    this.f1(p, 'BgPressure', K_BACKGROUND_PRESSURE);
    this.f2(p, 'PelletCentre', m.ballX, m.ballY);
    this.f1(p, 'PelletDensity', K_PELLET_DENSITY);
    this.f1(p, 'PelletRadius', K_PELLET_RADIUS);
    this.i1(p, 'CarryCoils', this.carryCoils ? 1 : 0);
    this.i1(p, 'PrevCoilCount', this.previousCoils.count);
    p.setArray('PrevCoilData', coilData(this.previousCoils), 3);
    this.f1(p, 'PrevGuideBz', this.previousCoils.guide);
    const s = this.state[this.current];
    this.withTextures([...s.textures, input], () => this.drawInto(this.state[1 - this.current]));
    this.current = 1 - this.current;
  }

  //-------------------------------------------------------------------------
  // `Clock()`: the reduction chain, then the clock.
  //-------------------------------------------------------------------------
  clockPass(substep, substeps, target) {
    const gl = this.gl;
    {
      const p = this.programs.reduce.use();
      this.samplers(p, ['Source']);
      let source = this.state[this.current].textures[2];
      let sw = this.grid.nx;
      let sh = this.grid.ny;
      let fromState = true;
      this.i1(p, 'Block', K_REDUCE_BLOCK);
      for (const level of this.reduction) {
        this.i2(p, 'SourceSize', sw, sh);
        this.i1(p, 'FromState', fromState ? 1 : 0);
        this.withTextures([source], () => this.drawInto(level));
        source = level.texture;
        sw = level.width;
        sh = level.height;
        fromState = false;
      }
    }

    const m = this.currentModel();
    const p = this.programs.clock.use();
    const next = 1 - this.clockIndex;
    this.samplers(p, ['Reduced', 'PrevClock']);
    const last = this.reduction[this.reduction.length - 1];
    this.i2(p, 'ReducedSize', last.width, last.height);
    this.i1(p, 'Substep', substep);
    this.i1(p, 'Substeps', substeps);
    this.f1(p, 'Target', target);
    this.f1(p, 'Courant', K_COURANT);
    this.f1(p, 'Dx', this.grid.dx);
    this.f1(p, 'Eta', m.eta);
    this.withTextures([last.texture, this.clock[this.clockIndex].texture], () => {
      gl.bindFramebuffer(gl.FRAMEBUFFER, this.clock[next].fbo);
      gl.viewport(0, 0, 2, 1);
      this.quad.draw();
    });
    this.clockIndex = next;
  }

  //-------------------------------------------------------------------------
  // `Substep()`: predict, the two face passes, the update.
  //-------------------------------------------------------------------------
  substep() {
    const m = this.currentModel();
    const s = this.state[this.current];
    const clk = this.clock[this.clockIndex].texture;
    {
      const p = this.programs.predict.use();
      this.setStateUniforms(p);
      this.samplers(p, ['StateA', 'StateB', 'StateC', 'StateD', 'ClockTexture']);
      this.withTextures([...s.textures, clk], () => this.drawInto(this.star));
    }
    {
      const p = this.programs.flux.use();
      this.setStateUniforms(p);
      this.samplers(p, ['StateA', 'StateB', 'StateC', 'StateD', 'StarA', 'StarB', 'StarC', 'StarD', 'ClockTexture']);
      this.f1(p, 'Eta', m.eta);
      this.withTextures([...s.textures, ...this.star.textures, clk], () => {
        this.i1(p, 'Direction', 0);
        this.drawInto(this.fluxX);
        this.i1(p, 'Direction', 1);
        this.drawInto(this.fluxY);
      });
    }
    {
      const p = this.programs.update.use();
      this.setStateUniforms(p);
      this.samplers(p, ['StateA', 'StateB', 'StateC', 'StateD', 'StarA', 'StarB', 'ClockTexture', 'FluxXA', 'FluxXB',
        'FluxXC', 'FluxXD', 'FluxYA', 'FluxYB', 'FluxYC', 'FluxYD']);
      this.f1(p, 'CoolingRate', m.cooling);
      this.f1(p, 'CoolingFloorT', K_BACKGROUND_PRESSURE / K_BACKGROUND_DENSITY);
      this.f1(p, 'GLMAlpha', K_GLM_ALPHA);
      this.i1(p, 'UseEntropy', 1);
      this.f1(p, 'EntropySwitch', K_ENTROPY_SWITCH);
      this.f1(p, 'SpongeEFolds', this.grid.ox > 0 ? K_SPONGE_EFOLDS : 0);
      this.f1(p, 'SpongePower', K_SPONGE_POWER);
      this.withTextures([...s.textures, this.star.textures[0], this.star.textures[1], clk,
        ...this.fluxX.textures, ...this.fluxY.textures], () => this.drawInto(this.state[1 - this.current]));
    }
    this.current = 1 - this.current;
  }

  /// `ReadClock()`: ( dt, done, c_h, floors ), and the entropy tally. A
  /// synchronous read, as glGetTexImage is in the plugin.
  readClock() {
    const gl = this.gl;
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.clock[this.clockIndex].fbo);
    gl.readPixels(0, 0, 2, 1, gl.RGBA, gl.FLOAT, this.readback);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this.lastEntropyCells = this.readback[4];
    return [this.readback[0], this.readback[1], this.readback[2], this.readback[3]];
  }

  //-------------------------------------------------------------------------
  // `Advance()`: plan the substeps from last frame's speed (or measure it now
  // after an Ignite or a pellet), then run them. dt itself lives on the GPU.
  //-------------------------------------------------------------------------
  advance(target, synchronous) {
    if (target <= 0) return;
    if (synchronous || !this.speedKnown) {
      this.clockPass(0, 1, 0.0);
      const values = this.readClock();
      this.lastSpeed = f32(2.0 * values[2]);
      this.speedKnown = true;
    }
    const stable = K_COURANT * this.grid.dx / Math.max(this.lastSpeed, 1e-12);
    let step = stable;
    const m = this.currentModel();
    if (m.eta > 0) step = Math.min(step, 0.2 * this.grid.dx * this.grid.dx / m.eta);
    const wanted = Math.ceil(1.1 * target / step);
    this.frameSubsteps = clamp(wanted, 1, K_MAX_SUBSTEPS);

    for (let k = 0; k < this.frameSubsteps; k += 1) {
      this.clockPass(k, this.frameSubsteps, target);
      this.substep();
    }
    this.lastTarget = target;
    this.clockPending = true;
  }

  //-------------------------------------------------------------------------
  // The light: `Emission()`, `Glow()`, `Potential()`.
  //-------------------------------------------------------------------------
  emissionPass() {
    const m = this.currentModel();
    const p = this.programs.emission.use();
    this.setStateUniforms(p);
    this.samplers(p, ['StateA', 'StateB', 'StateC', 'StateD']);
    this.i1(p, 'View', this.option('view', 7));
    this.f1(p, 'Tint', clamp01(this.P('tint')));
    this.i1(p, 'RampKind', this.option('ramp', 2));
    this.f1(p, 'TempReference', m.pressure);
    this.f1(p, 'FieldReference', m.field);
    this.withTextures(this.state[this.current].textures, () => this.drawInto(this.emission));
  }

  glowPass() {
    this.glowFraction = glowFromParam(this.P('glow'));
    if (this.glowFraction <= 0) return;

    let source = this.emission.texture;
    {
      const p = this.programs.downsample.use();
      this.samplers(p, ['Source']);
      for (const level of this.glowReduce) {
        this.withTextures([source], () => this.drawInto(level));
        source = level.texture;
      }
    }

    const p = this.programs.blur.use();
    this.samplers(p, ['Source']);
    this.i2(p, 'Size', this.glowWidth, this.glowHeight);
    let previous = 0;
    for (let stage = 0; stage < K_GLOW_STAGES; stage += 1) {
      const sigma = K_GLOW_SIGMA[stage] * this.glowHeight;
      const step = Math.sqrt(Math.max(sigma * sigma - previous * previous, 1e-6));
      previous = sigma;
      const [taps, weights] = gaussianWeights(step);
      this.i1(p, 'Taps', taps);
      p.setArray('Weights', weights, 1);
      const passes = [[source, this.glowTemp, 1, 0], [this.glowTemp.texture, this.glow[stage], 0, 1]];
      for (const [from, into, dx, dy] of passes) {
        this.i2(p, 'Direction', dx, dy);
        this.withTextures([from], () => this.drawInto(into));
      }
      source = this.glow[stage].texture;
    }
  }

  potentialPass() {
    const residual = this.programs.residual;
    const smooth = this.programs.smooth;
    const prolong = this.programs.prolong;
    const data = coilData(this.coils);
    const stateB = this.state[this.current].textures[1];

    const common = (p, level, finest) => {
      this.i2(p, 'Size', level.nx, level.ny);
      this.f1(p, 'Spacing', level.spacing);
      this.i1(p, 'VacuumGhosts', finest ? 1 : 0);
      this.i1(p, 'CoilCount', this.coils.count);
      p.setArray('CoilData', data, 3);
      this.samplers(p, ['Solution', 'Rhs', 'StateB', 'Coarse']);
    };
    const relax = (level, finest, sweeps) => {
      smooth.use();
      common(smooth, level, finest);
      for (let i = 0; i < sweeps; i += 1) {
        this.withTextures([level.solution[level.current].texture, level.rhs.texture],
          () => this.drawInto(level.solution[1 - level.current]));
        level.current = 1 - level.current;
      }
    };

    const top = this.levels[0];
    residual.use();
    common(residual, top, true);
    this.i1(residual, 'Mode', 0);
    this.withTextures([top.residual.texture, top.residual.texture, stateB, top.residual.texture],
      () => this.drawInto(top.rhs));
    if (!this.havePotential) {
      top.solution[0].clearTo(0, 0, 0, 0);
      top.solution[1].clearTo(0, 0, 0, 0);
      relax(top, true, 20);
      this.havePotential = true;
    }

    const count = this.levels.length;
    for (let l = 0; l < count; l += 1) {
      const level = this.levels[l];
      if (l > 0) {
        level.solution[0].clearTo(0, 0, 0, 0);
        level.solution[1].clearTo(0, 0, 0, 0);
        level.current = 0;
      }
      const coarsest = l === count - 1;
      relax(level, l === 0, coarsest ? 40 : 2);
      if (coarsest) break;

      residual.use();
      common(residual, level, l === 0);
      this.i1(residual, 'Mode', 1);
      this.withTextures([level.solution[level.current].texture, level.rhs.texture, stateB],
        () => this.drawInto(level.residual));
      const down = this.programs.downsample.use();
      this.samplers(down, ['Source']);
      this.withTextures([level.residual.texture], () => this.drawInto(this.levels[l + 1].rhs));
    }
    for (let l = count - 2; l >= 0; l -= 1) {
      const level = this.levels[l];
      prolong.use();
      common(prolong, level, l === 0);
      const coarse = this.levels[l + 1];
      this.withTextures([level.solution[level.current].texture, level.rhs.texture, stateB,
        coarse.solution[coarse.current].texture], () => this.drawInto(level.solution[1 - level.current]));
      level.current = 1 - level.current;
      relax(level, l === 0, 2);
    }
  }

  //-------------------------------------------------------------------------
  // `ProcessOpenGL()`.
  //-------------------------------------------------------------------------
  render({ input, params, width, height, time }) {
    const gl = this.gl;
    this.params = params;
    gl.disable(gl.BLEND);
    this.takeEvents(params);

    // Time. This page's clock is already in seconds, so the plugin's vote on
    // the host's unit has nothing to decide. Restart sends the clock back to
    // zero, which the page reads as "from the top": a fresh ball.
    if (this.lastNow >= 0 && time < this.lastNow) {
      this.ignitePending = true;
      this.lastNow = -1;
    }
    const hostDt = this.lastNow >= 0 ? clamp(time - this.lastNow, 0, K_MAX_FRAME_DELTA) : 0;
    this.lastNow = time;

    // Last frame's clock, read now rather than then.
    if (this.clockPending) {
      const values = this.readClock();
      this.lastSpeed = f32(2.0 * values[2]);
      telemetry.floors += values[3];
      telemetry.entropy = this.lastEntropyCells;
      telemetry.done = values[1];
      telemetry.target = this.lastTarget;
      this.clockPending = false;
      if (values[1] < 0.999 * this.lastTarget) telemetry.capped += 1;
    }

    // The grid, and every allocation.
    const detail = this.option('detail', DETAIL_CELLS.length);
    const boundaryNow = this.option('boundary', 2);
    const margin = boundaryNow === BOUNDARY_OPEN ? Math.round(K_MARGIN_FRACTION * DETAIL_CELLS[detail]) : 0;
    const wanted = chooseGrid(input.width, input.height, DETAIL_CELLS[detail], margin);
    this.ensureBuffers(wanted);

    // The coils, as of this frame.
    const m = this.currentModel();
    const simDt = hostDt * speedFromParam(this.P('speed'));
    this.spinAngle += m.spin * simDt;
    this.spinAngle -= Math.round(this.spinAngle / (2.0 * K_PI)) * 2.0 * K_PI;
    const wantStrength = this.P('quench') >= 0.5 ? 0.0 : 1.0;
    this.coilStrength += (wantStrength - this.coilStrength) * (1.0 - Math.exp(-simDt / K_QUENCH_TIME));
    this.previousCoils = this.coils;
    this.coils = makeCoils(m.poles, m.coilRadius, this.spinAngle, m.field, m.guide, this.coilStrength,
      0.5 * this.grid.lx, 0.5 * this.grid.ly);
    this.carryCoils = (m.boundary === BOUNDARY_WALL || m.boundary === BOUNDARY_OPEN)
      && this.previousCoils.count === this.coils.count && !this.ignitePending && this.grid.nx > 0;
    this.drive.advance(simDt, m.driveScale, m.drive);

    // Events, then the plasma.
    let synchronous = false;
    if (this.ignitePending) {
      this.ignite(input.texture);
      this.ignitePending = false;
      synchronous = true;
    }
    const pellets = this.pelletPresses;
    this.pelletPresses = 0;
    this.sources(input.texture, hostDt, simDt, pellets);
    if (pellets > 0) synchronous = true;

    this.advance(simDt, synchronous);
    this.simTime += simDt;

    // The light.
    this.emissionPass();
    const view = this.option('view', 7);
    if (view === 0) this.glowPass();
    const lines = clamp01(this.P('fieldLines'));
    let lineSpacing = 0;
    if (view === 0 && lines > 0 && this.coils.count > 0) {
      this.potentialPass();
      const cx = 0.5 * this.grid.lx;
      const cy = 0.5 * this.grid.ly;
      const centre = vacuumPotential(this.coils, cx, cy);
      let span = 0;
      for (let k = 0; k < 360; k += 1) {
        const a = 2.0 * K_PI * k / 360.0;
        span = Math.max(span, Math.abs(vacuumPotential(this.coils, cx + 0.5 * Math.cos(a), cy + 0.5 * Math.sin(a)) - centre));
      }
      lineSpacing = span / lineCountFromParam(this.P('lineCount'));
    }

    // Composite, onto the canvas.
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, width, height);
    const p = this.programs.composite.use();
    this.samplers(p, ['InputTexture', 'Emission', 'Glow0', 'Glow1', 'Glow2', 'Potential']);
    this.f2(p, 'MaxUV', 1, 1);
    const exact = input.width === width && input.height === height;
    this.i1(p, 'ExactInput', exact ? 1 : 0);
    this.i1(p, 'View', view);
    const ev = exposureFromParam(this.P('exposure'));
    this.f1(p, 'Gain', Math.pow(2.0, ev) / K_EMISSION_REFERENCE);
    const glowAmount = view === 0 ? this.glowFraction : 0;
    this.f1(p, 'GlowAmount', glowAmount);
    this.f1(p, 'CoreWeight', 1.0 - glowAmount);
    p.set('GlowShare', K_GLOW_SHARE[0], K_GLOW_SHARE[1], K_GLOW_SHARE[2]);
    this.f1(p, 'LineAmount', lines);
    this.f1(p, 'LineSpacing', lineSpacing);
    const g = this.grid;
    p.set('PotentialMap', g.fx / g.nx, g.fy / g.ny, g.ox / g.nx, g.oy / g.ny);
    this.f1(p, 'MixAmount', clamp01(this.P('mix')));
    const potential = this.havePotential ? this.levels[0].solution[this.levels[0].current].texture : this.emission.texture;
    this.withTextures([input.texture, this.emission.texture, this.glow[0].texture, this.glow[1].texture,
      this.glow[2].texture, potential], () => this.quad.draw());

    telemetry.grid = `${g.fx}×${g.fy}${g.ox ? ` (+${g.ox} margin)` : ''}`;
    telemetry.cells = g.nx * g.ny;
    telemetry.substeps = this.frameSubsteps;
    telemetry.frames += 1;
    telemetry.speed = speedFromParam(this.P('speed'));
  }
}

function coilData(coils) {
  const data = new Float32Array(K_MAX_COILS * 3);
  for (let k = 0; k < coils.count; k += 1) {
    data[k * 3 + 0] = coils.x[k];
    data[k * 3 + 1] = coils.y[k];
    data[k * 3 + 2] = coils.current[k];
  }
  return data;
}

/// `gaussianWeights()`: taps, and the weights for |offset| = 0..taps,
/// normalised over the full, symmetric kernel.
function gaussianWeights(sigma) {
  const taps = clamp(Math.ceil(3.0 * sigma), 1, K_MAX_BLUR_TAPS);
  const weights = new Float32Array(K_MAX_BLUR_TAPS + 1);
  const w = [];
  let sum = 0;
  for (let k = 0; k <= taps; k += 1) {
    w.push(Math.exp(-0.5 * k * k / Math.max(sigma * sigma, 1e-12)));
    sum += k === 0 ? w[k] : 2.0 * w[k];
  }
  for (let k = 0; k <= taps; k += 1) weights[k] = w[k] / sum;
  return [taps, weights];
}

//===========================================================================
// The parameters, from the constructor in Containment.cpp, in its order and
// groups. Audio, Audio Heat and Audio Pellets are left out, and so is the
// About block: see the note at the top.
//===========================================================================

const std = (id, name, def, group, extra = {}) => ({ id, name, type: 'standard', default: def, group, ...extra });
const button = (id, name, group, hint) => ({ id, name, type: 'boolean', default: 0, group, hint });
const opt = (id, name, def, group, elements, hint) => ({ id, name, type: 'option', default: def, group, elements, hint });
const fixed = (v, d) => v.toFixed(d);

const PARAMS = [
  opt('preset', 'Preset', 0, 'Preset', ['Custom', ...PRESETS.map((p) => p.name)],
    'The plugin’s own presets. Custom means the controls are the truth; any other row is laid over them while it is chosen (the sliders keep showing your values, as Resolume’s do), and choosing one re-ignites. Clip Orb is the constructor’s defaults.'),

  button('ignite', 'Ignite', 'Ball', 'Lays down a fresh ball. An event button in the plugin; here a toggle the renderer releases on the frame it acts, which is why it blinks.'),
  std('ballSize', 'Ball Size', 0.389, 'Ball', {
    display: (v) => `${fixed(ballSizeFromParam(v), 3)} frame heights`,
    hint: 'The ball’s radius: the 1/e radius of the Gaussian, or the top hat’s edge. Ignite lays it down at this size; Fuel tops up the same footprint.' }),
  std('ballX', 'Ball X', 0.5, 'Ball', { display: (v) => `${(v * 100).toFixed(0)}% across` }),
  std('ballY', 'Ball Y', 0.5, 'Ball', { display: (v) => `${(v * 100).toFixed(0)}% up` }),
  std('temperature', 'Temperature', 0.315, 'Ball', {
    display: (v) => `β₀ ${fixed(temperatureFromParam(v), 2)}`,
    hint: 'The ball’s β against the reference field, 0.1 to 50 geometrically: p₀ = β₀/2. Field confines the same ball harder rather than rescaling it.' }),
  opt('profile', 'Profile', 0, 'Ball', ['Gaussian', 'Top Hat']),
  button('pellet', 'Pellet', 'Ball', 'Drops a cold, dense pellet into the ball’s centre. An event button in the plugin; here a toggle the renderer releases on the frame it acts.'),
  std('feed', 'Feed', 0.1, 'Ball', {
    display: (v) => `${(feedFromParam(v) * 100).toFixed(0)}% a frame`,
    hint: 'How far the plasma’s colour moves towards the clip each frame (at 60 fps), inside the ball.' }),
  std('clipHeats', 'Clip Heats', 0.0, 'Ball', {
    display: (v) => fixed(clipHeatsFromParam(v), 2),
    hint: 'At 1 the brightest pixel starts at twice the ball’s pressure and the darkest at none, and Feed keeps heating it by the clip’s light.' }),
  std('fuel', 'Fuel', 0.15, 'Ball', {
    display: (v) => `${fixed(fuelFromParam(v), 2)} per τ_A`,
    hint: 'A steady gas puff that tops the ball’s footprint back up towards its ignition profile. It only adds, so it holds a leaking ball up without pinning it.' }),

  std('field', 'Field', 0.25, 'Bottle', {
    display: (v) => `${fixed(fieldFromParam(v), 2)} B_ref`,
    hint: 'The bottle’s strength: |B| at radius 0.5 in the first gap between coils.' }),
  std('guideField', 'Guide Field', 0.7, 'Bottle', {
    display: (v) => `${fixed(guideFieldFromParam(v), 2)} × Field`,
    hint: 'The uniform axial field Bz. It holds the ball; the cusp shapes it.' }),
  opt('poles', 'Poles', 2, 'Bottle', ['0', '4', '6', '8', '12'], 'Line currents of alternating sign round the frame. 0 is the guide field alone.'),
  std('coilRadius', 'Coil Radius', 0.105, 'Bottle', {
    display: (v) => `${fixed(coilRadiusFromParam(v), 2)} frame heights`,
    hint: 'Never inside the frame: the plugin pushes the coils out to 1.05 × the grid’s half-diagonal.' }),
  std('coilSpin', 'Coil Spin', 0.55, 'Bottle', {
    display: (v) => `${fixed(coilSpinFromParam(v), 2)} rad per τ_A`,
    hint: 'The coils turn; the change is carried through the vessel.' }),
  std('curvature', 'Curvature', 0.316, 'Bottle', {
    display: (v) => `g_eff ${fixed(curvatureFromParam(v), 2)}`,
    hint: 'An effective gravity, radial from the ball, standing in for bad curvature. It drives interchange (magnetic Rayleigh–Taylor) at the ball’s edge.' }),
  button('quench', 'Quench', 'Bottle', 'The coils lose their current (τ 0.15 τ_A) and the ball free-expands into a fireball. They come back when it is released.'),
  opt('boundary', 'Boundary', 1, 'Bottle', ['Open', 'Wall'],
    'Wall: a conducting, no-slip vessel. Open: a window onto a bigger bottle, with an absorbing margin simulated round the frame and never shown — a third more cells.'),

  std('speed', 'Speed', paramFromSpeed(0.3), 'Plasma', {
    display: (v) => (speedFromParam(v) === 0 ? 'frozen' : `${fixed(speedFromParam(v), 3)} τ_A per second`),
    hint: 'Alfvén crossing times per second, 0.02 to 2, and exactly zero at the bottom.' }),
  std('resistivity', 'Resistivity', 0.0, 'Plasma', {
    display: (v) => (resistivityFromParam(v) === 0 ? 'ideal' : `η ${resistivityFromParam(v).toExponential(1)}`) }),
  std('cooling', 'Cooling', 0.0, 'Plasma', {
    display: (v) => `Λ ${fixed(coolingFromParam(v), 2)}`,
    hint: 'Bremsstrahlung cooling, dE/dt = −Λ ρ² √T.' }),
  opt('detail', 'Detail', DEMO_DETAIL, 'Plasma', ['128', '256', '512', '1024'],
    'Grid cells on the frame’s short side. The plugin’s default is 256; this page starts at 128 so a browser keeps up. Changing it re-ignites.'),
  std('drive', 'Drive', 0.274, 'Plasma', {
    display: (v) => `${fixed(driveFromParam(v), 2)} rms`,
    hint: 'The stirring: twelve Fourier modes of a stream function with wandering (Ornstein–Uhlenbeck) amplitudes, windowed about the ball. Divergence-free.' }),
  std('driveScale', 'Drive Scale', 0.442, 'Plasma', {
    display: (v) => `${fixed(driveScaleFromParam(v), 3)} frame heights`,
    hint: 'The stirring’s wavelength.' }),

  std('exposure', 'Exposure', 0.4167, 'Light', {
    display: (v) => `${exposureFromParam(v) >= 0 ? '+' : ''}${fixed(exposureFromParam(v), 1)} EV` }),
  std('tint', 'Temperature Tint', 0.0, 'Light', {
    display: (v) => fixed(clamp01(v), 2),
    hint: '0 keeps the clip’s own colours; 1 colours the plasma by its temperature on the Ramp.' }),
  opt('ramp', 'Ramp', 0, 'Light', ['Hot', 'Aurora'], 'Hot: deep red to violet to white. Aurora: the 557.7 nm oxygen green, going white-green when hot.'),
  std('glow', 'Glow', 0.333, 'Light', {
    display: (v) => `${(glowFromParam(v) * 100).toFixed(0)}% into the glare`,
    hint: 'The fraction of the light the camera’s glare moves out of the core. It moves light; it never adds any.' }),
  std('fieldLines', 'Field Lines', 0.0, 'Light', {
    display: (v) => fixed(clamp01(v), 2),
    hint: 'Contours of A_z, solved every frame by a multigrid V-cycle, lit by the plasma on them.' }),
  std('lineCount', 'Field Line Count', 0.273, 'Light', { display: (v) => `${lineCountFromParam(v)} lines` }),
  opt('view', 'View', 0, 'Light', ['Picture', 'Density', 'Pressure', '|B|', 'Beta', 'Speed', 'Div B'],
    'Picture is the effect; the rest show one field of the state.'),
  std('mix', 'Mix', 1.0, 'Light', { display: (v) => `${(v * 100).toFixed(0)}%` }),
];

const demo = mountDemo({
  name: 'Containment',
  pluginId: 'CT01',
  tagline: 'A ball of plasma in a magnetic bottle. The clip is what the plasma is made of, and the plasma obeys 2.5-D compressible ideal MHD — MUSCL–Hancock, HLLD, GLM divergence cleaning — so its expansion, its ringing, its interchange fingers and its leaks through the cusps are the equations, not an animation.',
  repo: 'https://github.com/stoatworks-labs/containment',
  blurb: 'It is Containment’s own GLSL solver, ported from the repository to WebGL2 and run for real on a generated clip in this page: the same passes, the same parameters, the same maths, with the plugin’s CPU half ported to JavaScript. It starts at Detail 128, not the plugin’s 256, so a browser keeps up.',

  // Four RGBA32F state textures, written four at a time; float read-back of
  // the clock. OES_texture_float_linear is required in createRenderer.
  needFloat: true,

  // The clip is what the plasma is made of: its colour, carried with the
  // mass. The bars make a vivid ball whose stripes show the flow; the scene
  // is darker, and the ball with it.
  sources: ['bars', 'scene', 'spot', 'grid', 'ramp', 'detail'],

  params: PARAMS,

  differences: [
    'Every piece of the plugin’s GLSL is here unedited, assembled into the same fourteen programs, and demo/tools/check_shaders.py proves both. Everything the CPU does is a hand port to JavaScript that nothing checks but a reader: the control conversions, the coils’ closed-form field and its normalisation, the grid choice, the PCG generator and the Ornstein–Uhlenbeck stirring, the presets and their override, the substep plan, the events, the Gaussian weights of the glow and the multigrid schedule.',
    'The grid. The plugin’s default is Detail 256 (456 × 256 cells on a 16:9 frame, about 22 substeps of seven dependent GPU passes a frame). This page defaults to Detail 128 (224 × 128) so a browser keeps up; 256, 512 and 1024 are on the dropdown and run the same solver. At 128 the plasma’s edge and fingers are coarser than the plugin draws them by default. Measured headless in Chrome on an Apple M4 Max: 128 holds 60 frames a second at 12–14 substeps, 256 about 50–60 at 26–31; a laptop or phone GPU will manage far less at 256.',
    'Substeps. The plugin plans each frame’s substeps from the last frame’s speed, read back one frame late, plus 10%, capped at 48; past the cap simulated time runs slow rather than unstable. The page does exactly that, including the synchronous read-back after an Ignite or a pellet, and the stats line under the picture says when the cap bites. A browser that cannot hold the frame rate simply runs fewer frames: each frame is still clamped to a tenth of a second, as in the plugin, so a slow machine sees the plasma evolve slowly, never differently.',
    'The audio side is not here. The plugin reads Resolume’s FFT buffer: Audio Heat heats the ball with the level and Audio Pellets drops a pellet on each onset. A browser has no such buffer, and asking for a microphone to demonstrate a video effect is not a trade worth making, so the three controls are absent rather than present and dead. The removal is exact: with no spectrum the plugin’s analyser reports no level and never fires.',
    'Ignite and Pellet are FF_TYPE_EVENT buttons in the plugin; the kit has no event type, so each is a toggle the renderer takes and releases on the frame it acts, which is what a host does with an event anyway. Quench is a plain toggle, as in the plugin.',
    'The Preset dropdown is the plugin’s, with its override model: a row is laid over the controls while it is chosen and the sliders keep showing your values, exactly as they do in Resolume. The kit’s own presets menu is not used.',
    'Float render targets. The solver needs EXT_color_buffer_float (the state is four 32-bit float textures, and the clock is read back as floats) and OES_texture_float_linear (the light and the potential are read with a linear filter). Both are core in the plugin’s GL 4.1 and optional in WebGL2; the page stops with a message if either is missing rather than render a plausible wrong plasma. The state is never held in half floats.',
    'The clock is this page’s, in seconds, so the plugin’s vote on whether the host sends seconds or milliseconds has nothing to decide. Restart lays down a fresh ball; the plugin has no such control (Ignite is its equivalent).',
    'The About block is not here; the links are in the header.',
  ],

  createRenderer: (gl, quad) => new ContainmentRenderer(gl, quad),
});

//===========================================================================
// The stats line. Not decoration: it is how a visitor tells a slow browser
// (fewer frames, the same plasma) from the substep cap biting (simulated time
// running slow, which the plugin logs through Diag). Skipped in embed mode.
//===========================================================================
if (demo && !new URLSearchParams(window.location.search).has('embed')) {
  const stage = document.querySelector('.stage');
  if (stage) {
    const line = document.createElement('p');
    line.className = 'stage__status';
    line.dataset.stats = '';
    stage.append(line);
    let lastFrames = 0;
    let lastAt = performance.now();
    setInterval(() => {
      const now = performance.now();
      const fps = ((telemetry.frames - lastFrames) * 1000) / Math.max(now - lastAt, 1);
      lastFrames = telemetry.frames;
      lastAt = now;
      if (!telemetry.grid) return;
      const slow = telemetry.target > 0 && telemetry.done < 0.999 * telemetry.target;
      const clock = telemetry.target > 0
        ? (slow
          ? `the substep cap bit: the last frame covered ${telemetry.done.toFixed(4)} of ${telemetry.target.toFixed(4)} τ_A, so simulated time is running slow`
          : `the last frame covered its ${telemetry.target.toFixed(4)} τ_A`)
        : 'the plasma is frozen';
      const branch = telemetry.cells && telemetry.substeps
        ? `, ${Math.round((100 * telemetry.entropy) / (telemetry.cells * telemetry.substeps))}% of cell-steps on the entropy branch`
        : '';
      line.textContent =
        `Grid ${telemetry.grid} cells, ${telemetry.substeps} substep${telemetry.substeps === 1 ? '' : 's'} last frame; ${clock}${branch}. `
        + `${fps.toFixed(0)} frames a second in this browser; floor clamps so far: ${telemetry.floors.toFixed(0)}.`;
    }, 500);
  }
}
