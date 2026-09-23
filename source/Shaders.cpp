#include "Shaders.h"

/*
    Every piece below is one raw string well under 15 KB (MSVC C2026 caps a
    single literal near 16 KB); the programs are assembled from them at run
    time by SourceFor(). Nothing here is a template: each piece is plain
    GLSL 4.10 that the assembled program includes once.
*/

namespace containment
{
namespace
{

const char* const kVersion = "#version 410 core\n";

//---------------------------------------------------------------------------
// The one vertex shader: the SDK's screen quad, straight through.
//---------------------------------------------------------------------------
const char* const kVertex = R"(
layout( location = 0 ) in vec4 vPosition;
layout( location = 1 ) in vec2 vUV;
out vec2 uv;
void main()
{
	gl_Position = vPosition;
	uv          = vUV;
}
)";

//---------------------------------------------------------------------------
// The grid, the gas and the coils: shared by every pass that touches state.
//---------------------------------------------------------------------------
const char* const kCommon = R"(
uniform ivec2 GridSize;   //cells
uniform float Dx;         //cell size, frame heights (square cells)
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
)";

//---------------------------------------------------------------------------
// Reconstruction, the Riemann solver and GLM.
//---------------------------------------------------------------------------
const char* const kRiemann = R"(
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
)";

//---------------------------------------------------------------------------
// The body forces, shared by the predictor and the update.
//---------------------------------------------------------------------------
const char* const kForces = R"(
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
)";

//---------------------------------------------------------------------------
// The step: the reduction, the clock, the predictor and the update.
//---------------------------------------------------------------------------
const char* const kReduceMain = R"(
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
)";

//The clock finishes the reduction itself (one fragment, a loop over the last
//level, at most a few hundred texels) -- every pass on this chain is a
//render pass the GPU must drain before the next, and at the default grid
//that latency, not the arithmetic, is what a substep costs.
const char* const kClockMain = R"(
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
)";

const char* const kPredictMain = R"(
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
)";

//A face on a conducting, reflecting wall: nothing crosses it. No mass, no
//picture, no energy, no flux of the field; the wall pushes back with the
//normal stress the ghost gave, which is all momentum does there.
const char* const kFluxMain = R"(
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
)";

const char* const kUpdateMain = R"(
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

	Q w  = toPrim( A, B, C, D );
	outA = A;
	outB = B;
	outC = vec4( C.x, signalSpeeds( w ).x, C.z, fired + 2.0 * onEntropy );
	outD = D;
}
)";

//---------------------------------------------------------------------------
// Laying down a ball, and the per-frame sources.
//---------------------------------------------------------------------------
const char* const kIgniteMain = R"(
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
	vec2 uv    = ( vec2( cell ) + 0.5 ) / vec2( GridSize );
	vec4 clip  = texture( InputTexture, uv * MaxUV );

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
)";

const char* const kSourcesMain = R"(
uniform sampler2D StateA;
uniform sampler2D StateB;
uniform sampler2D StateC;
uniform sampler2D StateD;
uniform sampler2D InputTexture;
uniform vec2 MaxUV;
uniform float FeedFraction; //this frame's step of the tracer towards the clip
uniform float ClipHeatRate; //dp for a white pixel inside the ball, this frame
uniform float AudioHeatRate;//dp inside the ball, this frame
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
	vec2 uv   = ( vec2( cell ) + 0.5 ) / vec2( GridSize );
	vec4 clip = clamp( texture( InputTexture, uv * MaxUV ), 0.0, 1.0 );
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
)";

//---------------------------------------------------------------------------
// The light, on the grid.
//---------------------------------------------------------------------------
const char* const kEmissionMain = R"(
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
	ivec2 cell = ivec2( gl_FragCoord.xy );
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
		//div B dx / |B|, the same measure `cttest --divb` asserts, times 10.
		vec3 bxp = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 1, 0 ) ).q1.yzw;
		vec3 bxm = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( -1, 0 ) ).q1.yzw;
		vec3 byp = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 0, 1 ) ).q1.yzw;
		vec3 bym = fetchPrim( StateA, StateB, StateC, StateD, cell + ivec2( 0, -1 ) ).q1.yzw;
		float div = 0.5 * ( ( bxp.x - bxm.x ) + ( byp.y - bym.y ) ) / max( length( b ), 1e-6 );
		col = diverging( clamp( 10.0 * div, -1.0, 1.0 ) );
	}
	fragColor = vec4( col, v );
}
)";

//Exact 2x box reduction: each output is the mean of its four inputs, so the
//sum over the grid, times four, is unchanged.
const char* const kDownsampleMain = R"(
uniform sampler2D Source;
out vec4 fragColor;
void main()
{
	ivec2 o = ivec2( gl_FragCoord.xy ) * 2;
	fragColor = 0.25 * ( texelFetch( Source, o, 0 ) + texelFetch( Source, o + ivec2( 1, 0 ), 0 )
	                     + texelFetch( Source, o + ivec2( 0, 1 ), 0 ) + texelFetch( Source, o + ivec2( 1, 1 ), 0 ) );
}
)";

//A normalised, discrete Gaussian with a half-sample mirror at the edges. A
//symmetric kernel reflected at the boundary is a symmetric operator, so its
//columns sum to one as its rows do: the blur moves light and never makes or
//loses any, right up to the frame's edge.
const char* const kBlurMain = R"(
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
)";

//---------------------------------------------------------------------------
// A_z by multigrid: del^2 A = -J_z, cell-centred, 5-point, weighted Jacobi.
//---------------------------------------------------------------------------
const char* const kPoissonCommon = R"(
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
)";

const char* const kResidualMain = R"(
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
)";

const char* const kSmoothMain = R"(
out vec4 fragColor;
void main()
{
	ivec2 c   = ivec2( gl_FragCoord.xy );
	float a   = texelFetch( Solution, c, 0 ).r;
	float jac = 0.25 * ( neighbours( c ) - Spacing * Spacing * texelFetch( Rhs, c, 0 ).r );
	fragColor = vec4( mix( a, jac, 0.8 ), 0.0, 0.0, 0.0 );
}
)";

const char* const kProlongMain = R"(
uniform sampler2D Coarse;
out vec4 fragColor;
void main()
{
	ivec2 c = ivec2( gl_FragCoord.xy );
	vec2 uv = ( vec2( c ) + 0.5 ) / vec2( Size );
	fragColor = vec4( texelFetch( Solution, c, 0 ).r + texture( Coarse, uv ).r, 0.0, 0.0, 0.0 );
}
)";

//---------------------------------------------------------------------------
// Onto the host's raster.
//---------------------------------------------------------------------------
const char* const kCompositeMain = R"(
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
			float f = texture( Potential, uv ).r / LineSpacing;
			float d = abs( fract( f + 0.5 ) - 0.5 ) / max( fwidth( f ), 1e-6 );
			light += LineAmount * 3.0 * exp( -d * d ) * e.rgb;
		}

		light *= Gain;
		col = vec3( shoulder( light.r ), shoulder( light.g ), shoulder( light.b ) );
	}
	else
		col = e.rgb;

	fragColor = mix( source, vec4( col, 1.0 ), MixAmount );
}
)";

std::string join( std::initializer_list< const char* > pieces )
{
	std::string out;
	for( const char* piece : pieces )
		out += piece;
	return out;
}

} // namespace

ProgramSource SourceFor( Program program )
{
	const std::string vertex = join( { kVersion, kVertex } );
	switch( program )
	{
	case Program::Reduce: return { "reduce", vertex, join( { kVersion, kReduceMain } ) };
	case Program::Clock: return { "clock", vertex, join( { kVersion, kClockMain } ) };
	case Program::Predict:
		return { "predict", vertex, join( { kVersion, kCommon, kRiemann, kForces, kPredictMain } ) };
	case Program::Flux:
		return { "flux", vertex, join( { kVersion, kCommon, kRiemann, kFluxMain } ) };
	case Program::Update:
		return { "update", vertex, join( { kVersion, kCommon, kRiemann, kForces, kUpdateMain } ) };
	case Program::Ignite: return { "ignite", vertex, join( { kVersion, kCommon, kIgniteMain } ) };
	case Program::Sources: return { "sources", vertex, join( { kVersion, kCommon, kSourcesMain } ) };
	case Program::Emission: return { "emission", vertex, join( { kVersion, kCommon, kEmissionMain } ) };
	case Program::Downsample: return { "downsample", vertex, join( { kVersion, kDownsampleMain } ) };
	case Program::Blur: return { "blur", vertex, join( { kVersion, kBlurMain } ) };
	case Program::Residual: return { "residual", vertex, join( { kVersion, kPoissonCommon, kResidualMain } ) };
	case Program::Smooth: return { "smooth", vertex, join( { kVersion, kPoissonCommon, kSmoothMain } ) };
	case Program::Prolong: return { "prolong", vertex, join( { kVersion, kPoissonCommon, kProlongMain } ) };
	case Program::Composite: return { "composite", vertex, join( { kVersion, kCompositeMain } ) };
	case Program::Count: break;
	}
	return { "none", "", "" };
}

std::vector< ProgramSource > AllSources()
{
	std::vector< ProgramSource > all;
	for( int i = 0; i < static_cast< int >( Program::Count ); ++i )
		all.push_back( SourceFor( static_cast< Program >( i ) ) );
	return all;
}

} // namespace containment
