#include "Checks.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <chrono>
#include <cstdio>
#include <functional>
#include <numeric>

namespace cttest
{
namespace
{
/// A physics rig: a square raster (the raster does not matter; the grid is
/// Detail cells on a side), the plasma quiet, and one frame rendered so the
/// grid exists. The frame's own Ignite is then overwritten by the check.
bool PhysicsRig( Rig& rig, int cells, const TestModel& model, int raster = 128 )
{
	rig.plugin.SetModelForTest( model );
	if( !rig.Init( raster, raster ) )
		return false;
	rig.Quiet();
	rig.Set( PT_DETAIL, DetailParam( cells ) );
	rig.Set( PT_SPEED, 0.0f );
	rig.Set( PT_GLOW, 0.0f );
	return rig.Render( 1 );
}

//===========================================================================
// A reference for Brio-Wu, independent of the plugin: 1-D ideal MHD in
// double precision, MUSCL-Hancock on primitive variables with a minmod
// limiter and the Rusanov (local Lax-Friedrichs) flux, at 8192 cells. It
// shares no code with the GPU scheme and a different Riemann solver, so
// agreement is not a scheme agreeing with itself.
//===========================================================================
struct Prim1
{
	double r, u, v, w, p, by, bz;
};

struct BrioWuReference
{
	int n;
	double gamma, bx;
	std::vector< Prim1 > cells;

	static void flux( const Prim1& s, double gamma, double bx, double f[ 7 ], double u[ 7 ] )
	{
		const double b2 = bx * bx + s.by * s.by + s.bz * s.bz;
		const double pt = s.p + 0.5 * b2;
		const double E  = s.p / ( gamma - 1.0 ) + 0.5 * s.r * ( s.u * s.u + s.v * s.v + s.w * s.w ) + 0.5 * b2;
		const double vb = s.u * bx + s.v * s.by + s.w * s.bz;
		u[ 0 ] = s.r;
		u[ 1 ] = s.r * s.u;
		u[ 2 ] = s.r * s.v;
		u[ 3 ] = s.r * s.w;
		u[ 4 ] = E;
		u[ 5 ] = s.by;
		u[ 6 ] = s.bz;
		f[ 0 ] = s.r * s.u;
		f[ 1 ] = s.r * s.u * s.u + pt - bx * bx;
		f[ 2 ] = s.r * s.u * s.v - bx * s.by;
		f[ 3 ] = s.r * s.u * s.w - bx * s.bz;
		f[ 4 ] = ( E + pt ) * s.u - bx * vb;
		f[ 5 ] = s.by * s.u - bx * s.v;
		f[ 6 ] = s.bz * s.u - bx * s.w;
	}

	static double fast( const Prim1& s, double gamma, double bx )
	{
		const double a2 = gamma * s.p / s.r;
		const double b2 = ( bx * bx + s.by * s.by + s.bz * s.bz ) / s.r;
		const double t2 = ( s.by * s.by + s.bz * s.bz ) / s.r;
		const double d  = a2 - b2;
		return std::sqrt( 0.5 * ( a2 + b2 + std::sqrt( d * d + 4.0 * a2 * t2 ) ) );
	}

	static Prim1 toPrim( const double u[ 7 ], double gamma, double bx )
	{
		Prim1 s;
		s.r  = u[ 0 ];
		s.u  = u[ 1 ] / s.r;
		s.v  = u[ 2 ] / s.r;
		s.w  = u[ 3 ] / s.r;
		s.by = u[ 5 ];
		s.bz = u[ 6 ];
		s.p  = ( gamma - 1.0 )
		      * ( u[ 4 ] - 0.5 * s.r * ( s.u * s.u + s.v * s.v + s.w * s.w )
		          - 0.5 * ( bx * bx + s.by * s.by + s.bz * s.bz ) );
		return s;
	}

	static double minmod( double a, double b )
	{
		return a * b <= 0.0 ? 0.0 : ( std::abs( a ) < std::abs( b ) ? a : b );
	}

	void Run( double tEnd )
	{
		const double dx = 1.0 / n;
		double t        = 0.0;
		std::vector< Prim1 > next( cells.size() );
		auto get = [ & ]( int i ) { return cells[ std::clamp( i, 0, n - 1 ) ]; };
		while( t < tEnd )
		{
			double smax = 0.0;
			for( const Prim1& s : cells )
				smax = std::max( smax, std::abs( s.u ) + fast( s, gamma, bx ) );
			double dt = 0.4 * dx / smax;
			if( t + dt > tEnd )
				dt = tEnd - t;

			//Slopes and Hancock half-step, in primitive-variable form via
			//the conserved face fluxes.
			std::vector< Prim1 > left( n ), right( n );
			for( int i = 0; i < n; ++i )
			{
				const Prim1 m = get( i - 1 ), c = get( i ), p = get( i + 1 );
				double cv[ 7 ] = { c.r, c.u, c.v, c.w, c.p, c.by, c.bz };
				double mv[ 7 ] = { m.r, m.u, m.v, m.w, m.p, m.by, m.bz };
				double pv[ 7 ] = { p.r, p.u, p.v, p.w, p.p, p.by, p.bz };
				double lo[ 7 ], hi[ 7 ];
				for( int k = 0; k < 7; ++k )
				{
					const double s = minmod( cv[ k ] - mv[ k ], pv[ k ] - cv[ k ] );
					lo[ k ]        = cv[ k ] - 0.5 * s;
					hi[ k ]        = cv[ k ] + 0.5 * s;
				}
				Prim1 L { lo[ 0 ], lo[ 1 ], lo[ 2 ], lo[ 3 ], lo[ 4 ], lo[ 5 ], lo[ 6 ] };
				Prim1 R { hi[ 0 ], hi[ 1 ], hi[ 2 ], hi[ 3 ], hi[ 4 ], hi[ 5 ], hi[ 6 ] };
				double fl[ 7 ], ul[ 7 ], fr[ 7 ], ur[ 7 ];
				flux( L, gamma, bx, fl, ul );
				flux( R, gamma, bx, fr, ur );
				for( int k = 0; k < 7; ++k )
				{
					ul[ k ] -= 0.5 * dt / dx * ( fr[ k ] - fl[ k ] );
					ur[ k ] -= 0.5 * dt / dx * ( fr[ k ] - fl[ k ] );
				}
				left[ i ]  = toPrim( ul, gamma, bx );
				right[ i ] = toPrim( ur, gamma, bx );
			}
			std::vector< std::array< double, 7 > > F( static_cast< size_t >( n + 1 ) );
			for( int f = 0; f <= n; ++f )
			{
				const Prim1 L = right[ std::clamp( f - 1, 0, n - 1 ) ];
				const Prim1 R = left[ std::clamp( f, 0, n - 1 ) ];
				double fl[ 7 ], ul[ 7 ], fr[ 7 ], ur[ 7 ];
				flux( L, gamma, bx, fl, ul );
				flux( R, gamma, bx, fr, ur );
				const double s = std::max( std::abs( L.u ) + fast( L, gamma, bx ), std::abs( R.u ) + fast( R, gamma, bx ) );
				for( int k = 0; k < 7; ++k )
					F[ f ][ k ] = 0.5 * ( fl[ k ] + fr[ k ] ) - 0.5 * s * ( ur[ k ] - ul[ k ] );
			}
			for( int i = 0; i < n; ++i )
			{
				double f[ 7 ], u[ 7 ];
				flux( cells[ i ], gamma, bx, f, u );
				for( int k = 0; k < 7; ++k )
					u[ k ] -= dt / dx * ( F[ i + 1 ][ k ] - F[ i ][ k ] );
				next[ i ] = toPrim( u, gamma, bx );
			}
			cells.swap( next );
			t += dt;
		}
	}
};

/// The features of a Brio-Wu profile at t = 0.1, located the same way in the
/// plugin's profile and the reference's.
struct Features
{
	double fastLeft;  ///< head of the left fast rarefaction
	double compound;  ///< where the tangential field changes sign
	double contact;   ///< steepest density gradient where pressure is continuous
	double slowShock; ///< steepest pressure gradient right of the contact
	double fastRight; ///< head of the right fast rarefaction
	double plateau1;  ///< density between the compound wave and the contact
	double plateau2;  ///< density between the contact and the slow shock
};

Features Locate( const std::vector< double >& x, const std::vector< double >& rho, const std::vector< double >& p,
                 const std::vector< double >& bt )
{
	const int n = static_cast< int >( x.size() );
	Features f {};
	//Heads: the first cell, from the outside in, that has moved off its
	//initial state by 1% of the rarefaction's own jump.
	for( int i = 0; i < n; ++i )
		if( std::abs( rho[ i ] - rho[ 0 ] ) > 1e-3 * rho[ 0 ] )
		{
			f.fastLeft = x[ i ];
			break;
		}
	for( int i = n - 1; i >= 0; --i )
		if( std::abs( rho[ i ] - rho[ n - 1 ] ) > 1e-3 * rho[ n - 1 ] )
		{
			f.fastRight = x[ i ];
			break;
		}
	//The compound wave: the tangential field's zero crossing, interpolated.
	for( int i = 0; i + 1 < n; ++i )
		if( bt[ i ] > 0.0 && bt[ i + 1 ] <= 0.0 && x[ i ] > 0.4 )
		{
			f.compound = x[ i ] + ( x[ i + 1 ] - x[ i ] ) * bt[ i ] / ( bt[ i ] - bt[ i + 1 ] );
			break;
		}
	//The contact and the slow shock: steepest gradients between the compound
	//wave and the right fast rarefaction. The contact is where density jumps
	//and pressure does not; the slow shock is where both do.
	double best = 0.0;
	for( int i = 1; i + 1 < n; ++i )
	{
		if( x[ i ] < f.compound + 0.01 || x[ i ] > 0.7 )
			continue;
		const double dr = std::abs( rho[ i + 1 ] - rho[ i - 1 ] );
		const double dp = std::abs( p[ i + 1 ] - p[ i - 1 ] );
		if( dr > best && dp < 0.25 * dr )
		{
			best      = dr;
			f.contact = x[ i ];
		}
	}
	best = 0.0;
	for( int i = 1; i + 1 < n; ++i )
	{
		if( x[ i ] < f.contact + 0.01 || x[ i ] > 0.75 )
			continue;
		const double dp = std::abs( p[ i + 1 ] - p[ i - 1 ] );
		if( dp > best )
		{
			best        = dp;
			f.slowShock = x[ i ];
		}
	}
	auto at = [ & ]( double where ) {
		int k = 0;
		while( k + 1 < n && x[ k + 1 ] < where )
			++k;
		return rho[ k ];
	};
	f.plateau1 = at( 0.5 * ( f.compound + f.contact ) );
	f.plateau2 = at( 0.5 * ( f.contact + f.slowShock ) );
	return f;
}

double FastSpeed( double r, double p, double bx, double by, double bz, double gamma )
{
	return BrioWuReference::fast( { r, 0, 0, 0, p, by, bz }, gamma, bx );
}
} // namespace

//===========================================================================
// --briowu
//===========================================================================
int RunBrioWu( const Perturb& perturb )
{
	Say( "\n=== briowu: the Brio-Wu shock tube (gamma 2, Bx 0.75), t = 0.1, along x and along y\n" );
	const double gamma = 2.0, bn = 0.75, tEnd = 0.1;

	//The reference, once.
	BrioWuReference ref;
	ref.n     = 8192;
	ref.gamma = gamma;
	ref.bx    = bn;
	for( int i = 0; i < ref.n; ++i )
	{
		const bool left = ( i + 0.5 ) / ref.n < 0.5;
		ref.cells.push_back( left ? Prim1 { 1.0, 0, 0, 0, 1.0, 1.0, 0 } : Prim1 { 0.125, 0, 0, 0, 0.1, -1.0, 0 } );
	}
	ref.Run( tEnd );
	std::vector< double > rx, rr, rp, rb;
	for( int i = 0; i < ref.n; ++i )
	{
		rx.push_back( ( i + 0.5 ) / ref.n );
		rr.push_back( ref.cells[ i ].r );
		rp.push_back( ref.cells[ i ].p );
		rb.push_back( ref.cells[ i ].by );
	}
	const Features R = Locate( rx, rr, rp, rb );

	//The two heads are known in closed form: the fast speed of each initial state.
	const double headL = 0.5 - tEnd * FastSpeed( 1.0, 1.0, bn, 1.0, 0.0, gamma );
	const double headR = 0.5 + tEnd * FastSpeed( 0.125, 0.1, bn, -1.0, 0.0, gamma );
	Say( "  reference (double, 8192 cells, Rusanov): fast %.4f  compound %.4f  contact %.4f  slow %.4f  fast %.4f\n",
	     R.fastLeft, R.compound, R.contact, R.slowShock, R.fastRight );
	Say( "  plateaus rho %.4f and %.4f; the fast heads in closed form at %.4f and %.4f\n", R.plateau1, R.plateau2,
	     headL, headR );
	Check( std::abs( R.fastLeft - headL ) < 4.0 / ref.n + 0.002 && std::abs( R.fastRight - headR ) < 4.0 / ref.n + 0.002,
	       fmt( "the reference's fast heads sit where the fast speeds put them (%.4f vs %.4f, %.4f vs %.4f)", R.fastLeft,
	            headL, R.fastRight, headR ) );

	TestModel model;
	model.gamma    = static_cast< float >( gamma );
	model.boundary = 3;//outflow, zero gradient in every variable
	model.solver   = perturb.hll ? 1 : 0;

	for( int cells : { 256, 512 } )
		for( int axis = 0; axis < 2; ++axis )
		{
			Rig rig;
			if( !PhysicsRig( rig, cells, model ) )
				return 1;
			const Grid& g = rig.plugin.CurrentGrid();
			StateBuilder sb( g.nx, g.ny, gamma );
			for( int j = 0; j < g.ny; ++j )
				for( int i = 0; i < g.nx; ++i )
				{
					const double s   = axis == 0 ? ( i + 0.5 ) / g.nx : ( j + 0.5 ) / g.ny;
					const bool left  = s < 0.5;
					const double rho = left ? 1.0 : 0.125, p = left ? 1.0 : 0.1, bt = left ? 1.0 : -1.0;
					if( axis == 0 )
						sb.Set( i, j, rho, 0, 0, 0, p, bn, bt, 0 );
					else
						sb.Set( i, j, rho, 0, 0, 0, p, bt, bn, 0 );
				}
			sb.Load( rig.plugin );
			const int steps = rig.plugin.StepForTest( tEnd );
			const Snapshot s = Snapshot::Take( rig.plugin, gamma );

			std::vector< double > x, r, p, b;
			const int n = axis == 0 ? g.nx : g.ny;
			for( int k = 0; k < n; ++k )
			{
				const int i = axis == 0 ? k : g.nx / 2;
				const int j = axis == 0 ? g.ny / 2 : k;
				x.push_back( ( k + 0.5 ) / n );
				r.push_back( s.Rho( i, j ) );
				p.push_back( s.P( i, j ) );
				b.push_back( axis == 0 ? s.By( i, j ) : s.Bx( i, j ) );
			}
			if( !s.Finite() || steps < 0 )
			{
				Check( false, fmt( "Detail %d along %c: the state went non-finite", cells, axis ? 'y' : 'x' ) );
				continue;
			}
			const Features F = Locate( x, r, p, b );
			//A discontinuity captured by a second-order scheme is spread over
			//two or three cells and its steepest point lies within one cell
			//of the true position; the heads, found by a 1% threshold, sit
			//up to three cells inside theirs. Three cells, derived from the
			//grid, is the bound for every feature.
			const double tol = 3.0 / n;
			const double worst = std::max( { std::abs( F.fastLeft - R.fastLeft ), std::abs( F.compound - R.compound ),
			                                 std::abs( F.contact - R.contact ), std::abs( F.slowShock - R.slowShock ),
			                                 std::abs( F.fastRight - R.fastRight ) } );
			Say( "  Detail %d along %c (%d substeps): fast %.4f  compound %.4f  contact %.4f  slow %.4f  fast %.4f\n",
			     cells, axis ? 'y' : 'x', steps, F.fastLeft, F.compound, F.contact, F.slowShock, F.fastRight );
			Check( worst <= tol, fmt( "Detail %d along %c: all five waves within 3 cells (%.4f) of the reference; worst %.4f "
			                          "(%.1f cells)",
			                          cells, axis ? 'y' : 'x', tol, worst, worst * n ) );
			//The plateaus: second-order schemes overshoot a slow shock by a
			//few percent over a cell or two, and the plateaus here are wide
			//(>= 20 cells at Detail 256), so their middles are clean. 2%.
			const double e1 = std::abs( F.plateau1 / R.plateau1 - 1.0 );
			const double e2 = std::abs( F.plateau2 / R.plateau2 - 1.0 );
			Check( e1 < 0.02 && e2 < 0.02, fmt( "Detail %d along %c: plateaus rho %.4f and %.4f against %.4f and %.4f "
			                                    "(%.2f%%, %.2f%%; bound 2%%)",
			                                    cells, axis ? 'y' : 'x', F.plateau1, F.plateau2, R.plateau1,
			                                    R.plateau2, 100 * e1, 100 * e2 ) );

			//The contact's width, 10% to 90% of its jump: HLLD resolves the
			//contact as a wave of its own; HLL has no contact wave and
			//smears it by its numerical diffusion. The bound is derived in
			//AGENTS.md: 4 cells at Detail 512.
			if( cells == 512 )
			{
				double lo = R.plateau1, hi = R.plateau2;
				const double a = std::min( lo, hi ) + 0.1 * std::abs( hi - lo );
				const double z = std::min( lo, hi ) + 0.9 * std::abs( hi - lo );
				int first = -1, last = -1;
				for( int k = 0; k < n; ++k )
				{
					if( std::abs( x[ k ] - F.contact ) > 0.05 )
						continue;
					if( r[ k ] <= z && r[ k ] >= a )
					{
						if( first < 0 )
							first = k;
						last = k;
					}
				}
				const int width = first < 0 ? 0 : last - first + 1;
				Check( width <= 4, fmt( "Detail 512 along %c: the contact is %d cells wide, 10%%-90%% (bound 4)",
				                        axis ? 'y' : 'x', width ) );
			}
		}
	return g_failures;
}


//===========================================================================
// --alfven
//===========================================================================
namespace
{
struct WaveFit
{
	double phaseBz, phasePerp, amplitude, correlation;
};

/// Project Bz and the in-plane perpendicular field onto cos and sin of k.x.
WaveFit FitWave( const Snapshot& s, double kx, double ky, double ex, double ey, double sqrtRho )
{
	double cz = 0, sz = 0, cp = 0, sp = 0, vb = 0, bb = 0, vv = 0;
	for( int j = 0; j < s.ny; ++j )
		for( int i = 0; i < s.nx; ++i )
		{
			const double x = ( i + 0.5 ) * s.dx, y = ( j + 0.5 ) * s.dx;
			const double ph = kx * x + ky * y;
			const double bz = s.Bz( i, j );
			const double bp = s.Bx( i, j ) * ex + s.By( i, j ) * ey;
			const double vp = s.U( i, j ) * ex + s.V( i, j ) * ey;
			cz += bz * std::cos( ph );
			sz += bz * std::sin( ph );
			cp += bp * std::cos( ph );
			sp += bp * std::sin( ph );
			vb += vp * bp;
			bb += bp * bp;
			vv += vp * vp;
		}
	const double n = static_cast< double >( s.nx ) * s.ny;
	WaveFit f;
	f.phaseBz     = std::atan2( sz, cz );
	f.phasePerp   = std::atan2( sp, cp );
	f.amplitude   = 2.0 * std::sqrt( cz * cz + sz * sz ) / n;
	f.correlation = vb / std::sqrt( std::max( bb * vv, 1e-300 ) );
	(void)sqrtRho;
	return f;
}
} // namespace

int RunAlfven( const Perturb& perturb )
{
	Say( "\n=== alfven: a circularly polarised Alfven wave on an oblique field, periodic box\n" );
	//Toth (2000)'s wave, on a unit box with k = 2 pi ( 1, 2 ): an exact
	//nonlinear solution whose only motion is the wave's. rho = 2 so that
	//B / sqrt( rho ) and B / rho differ by 41%.
	const double rho = 2.0, p = 0.1, b0 = 1.0, db = 0.1;
	const double kx = 2.0 * kPi, ky = 4.0 * kPi, k = std::sqrt( kx * kx + ky * ky );
	const double nx = kx / k, ny = ky / k;//along the wave
	const double ex = -ny, ey = nx;       //in-plane, across it
	const double vA = perturb.alfvenOverRho ? b0 / rho : b0 / std::sqrt( rho );
	const double trueVA = b0 / std::sqrt( rho );
	const double period = 2.0 * kPi / ( k * trueVA );

	TestModel model;
	model.boundary = 2;
	double loss[ 2 ] = {};
	int index = 0;
	for( int cells : { 128, 256 } )
	{
		Rig rig;
		if( !PhysicsRig( rig, cells, model ) )
			return 1;
		const Grid& g = rig.plugin.CurrentGrid();
		StateBuilder sb( g.nx, g.ny, kGamma );
		for( int j = 0; j < g.ny; ++j )
			for( int i = 0; i < g.nx; ++i )
			{
				//The cell average of a sinusoid is its centre value times
				//sinc( k dx / 2 ) per axis; loading averages keeps the
				//initial amplitude exact to second order and beyond.
				const double x = ( i + 0.5 ) * g.dx, y = ( j + 0.5 ) * g.dx;
				const double ph = kx * x + ky * y;
				const double avg = ( std::sin( 0.5 * kx * g.dx ) / ( 0.5 * kx * g.dx ) ) * ( std::sin( 0.5 * ky * g.dx ) / ( 0.5 * ky * g.dx ) );
				const double bp = db * std::sin( ph ) * avg, bz = db * std::cos( ph ) * avg;
				//A wave travelling along +n: delta v = -delta B / sqrt( rho ).
				const double vp = -bp / std::sqrt( rho ), vz = -bz / std::sqrt( rho );
				sb.Set( i, j, rho, vp * ex, vp * ey, vz, p, b0 * nx + bp * ex, b0 * ny + bp * ey, bz );
			}
		sb.Load( rig.plugin );
		const Snapshot s0 = Snapshot::Take( rig.plugin );
		const WaveFit f0  = FitWave( s0, kx, ky, ex, ey, std::sqrt( rho ) );

		//Four quarters of a period, so the phase can be unwrapped.
		double travelled = 0.0, last = f0.phaseBz;
		int steps        = 0;
		WaveFit f        = f0;
		for( int q = 0; q < 4; ++q )
		{
			steps += rig.plugin.StepForTest( 0.25 * period );
			const Snapshot s = Snapshot::Take( rig.plugin );
			f               = FitWave( s, kx, ky, ex, ey, std::sqrt( rho ) );
			double d        = f.phaseBz - last;
			while( d < -kPi )
				d += 2.0 * kPi;
			while( d > kPi )
				d -= 2.0 * kPi;
			travelled += d;
			last = f.phaseBz;
		}
		//Bz ~ cos( k.x - w t ): the fitted phase advances by w t.
		const double speed = travelled / ( k * period );
		//Second-order dispersion: the phase error of a centred second-order
		//scheme is O( ( k dx )^2 ); the bound is ( k dx )^2 itself.
		const double kdx   = k * g.dx;
		const double tol   = kdx * kdx;
		Check( std::abs( speed / vA - 1.0 ) < tol,
		       fmt( "Detail %d: phase speed %.5f against B/sqrt(rho) %.5f%s (error %.3f%%, bound (k dx)^2 = %.3f%%; %d substeps)",
		            cells, speed, trueVA, perturb.alfvenOverRho ? ", expected B/rho" : "",
		            100 * std::abs( speed / vA - 1.0 ), 100 * tol, steps ) );

		//Mode rotation: in a circularly polarised wave B_perp is a quarter
		//turn behind Bz (sin against cos) whatever the wave has travelled, and
		//v is anti-parallel to B_perp for a wave going along +B.
		double lead = f.phaseBz - f.phasePerp;
		while( lead < -kPi )
			lead += 2.0 * kPi;
		while( lead > kPi )
			lead -= 2.0 * kPi;
		Check( std::abs( lead + 0.5 * kPi ) < 0.02 && f.correlation < -0.99,
		       fmt( "Detail %d: B_perp a quarter turn from Bz (%.4f rad, -1.5708 expected), v against B_perp (corr %.4f)",
		            cells, lead, f.correlation ) );
		loss[ index++ ] = 1.0 - f.amplitude / f0.amplitude;
		Say( "  Detail %d: amplitude after one period %.5f of %.5f (lost %.3f%%)\n", cells, f.amplitude, f0.amplitude,
		     100 * loss[ index - 1 ] );
	}
	//Second order: the loss falls at least as 2^1.5 per doubling. (A smooth
	//wave is second order away from its extrema; the MC limiter clips each
	//extremum to first order over a few cells, so the global order of the
	//amplitude error lies between 1.5 and 2.)
	const double order = std::log( loss[ 0 ] / std::max( loss[ 1 ], 1e-12 ) ) / std::log( 2.0 );
	Check( order >= 1.5 && loss[ 1 ] > 0.0, fmt( "amplitude error converges at order %.2f (Detail 128 lost %.3f%%, 256 lost %.3f%%; bound >= 1.5)",
	                           order, 100 * loss[ 0 ], 100 * loss[ 1 ] ) );
	return g_failures;
}


//===========================================================================
// --conserve
//===========================================================================
namespace
{
struct Totals
{
	double mass, energy, absEnergy;
};
Totals Total( const Snapshot& s )
{
	Totals t { 0.0, 0.0, 0.0 };
	for( int j = 0; j < s.ny; ++j )
		for( int i = 0; i < s.nx; ++i )
		{
			t.mass += s.Rho( i, j );
			t.energy += s.Energy( i, j );
			t.absEnergy += std::abs( s.Energy( i, j ) );
		}
	return t;
}
} // namespace

int RunConserve( const Perturb& perturb )
{
	Say( "\n=== conserve: Wall, no cooling, no resistivity, no feed, no forcing: mass and energy over 600 frames\n" );
	//Two bottles. The guide field alone keeps the in-plane field zero, so GLM
	//has nothing to clean and beta stays well above the dual-energy switch:
	//there the scheme is conservative in every term and both totals must hold
	//to round-off. The six-pole cusp adds GLM's energy source and the switch
	//(both deliberately non-conservative, see AGENTS.md); there mass must
	//still hold to round-off and the energy drift is measured and reported.
	for( int poles : { 0, 6 } )
	{
		Rig rig;
		rig.Init( 256, 144 );
		rig.Quiet();
		rig.Set( PT_BOUNDARY, 1.0f );
		rig.Set( PT_POLES, PolesParam( poles ) );
		//A guide field of 0.3 keeps the background's beta (0.11) above the
		//switch's 0.02 everywhere, so in the guide-field bottle no cell may
		//take the entropy branch -- and none must, for the check to be of the
		//fluxes alone.
		rig.Set( PT_GUIDE_FIELD, GuideParam( 0.3 ) );
		rig.Set( PT_GLOW, 0.0f );
		if( perturb.conserveCooling )
			rig.Set( PT_COOLING, 0.3f );
		rig.Render( 1 );
		const Totals t0      = Total( Snapshot::Take( rig.plugin ) );
		const long long s0   = rig.plugin.SubstepsTaken();
		rig.Render( 600 );
		const Snapshot s     = Snapshot::Take( rig.plugin );
		const Totals t1      = Total( s );
		const long long steps = rig.plugin.SubstepsTaken() - s0;
		//Each substep rounds every cell's stored value a few times (three is
		//generous: the flux difference, the update, one source), each by at
		//most u = 2^-24 of that cell's own value. Summed over the cells that
		//is 3u of the total per step -- the cell count cancels, because the
		//error is relative to each cell and the check to their sum -- and the
		//worst case adds up linearly over the steps.
		const double u      = std::ldexp( 1.0, -24 );
		const double bound  = 3.0 * u * static_cast< double >( steps );
		const double dMass  = std::abs( t1.mass / t0.mass - 1.0 );
		const double dE     = std::abs( t1.energy - t0.energy ) / t0.absEnergy;
		Check( s.Finite() && dMass <= bound,
		       fmt( "%s: mass changed by %.2e of itself over %lld substeps (bound 3 u steps = %.2e)",
		            poles ? "six-pole cusp" : "guide field", dMass, steps, bound ) );
		const double switched = rig.plugin.EntropyCells() / ( static_cast< double >( s.nx ) * s.ny * steps );
		if( poles == 0 )
			Check( dE <= bound && rig.plugin.FloorHits() == 0.0 && rig.plugin.EntropyCells() == 0.0,
			       fmt( "guide field: energy changed by %.2e of itself (bound %.2e); floors %.0f, entropy-branch cells %.0f "
			            "(both must be 0)",
			            dE, bound, rig.plugin.FloorHits(), rig.plugin.EntropyCells() ) );
		else
			Say( "  six-pole cusp: energy changed by %.2e of itself -- GLM's energy source and the dual-energy switch "
			     "(on %.2f%% of cell-steps), not the fluxes (reported, not asserted)\n", dE, 100.0 * switched );
	}
	return g_failures;
}

//===========================================================================
// --divb
//===========================================================================
namespace
{
/// max and rms of |div B| dx / |B| over the grid, central differences, the
/// field's magnitude guarded below by a tenth of its rms (a cusp's null has
/// |B| = 0, where the ratio means nothing).
void DivergenceOf( const Snapshot& s, double& maxRatio, double& rmsRatio )
{
	double b2 = 0.0;
	for( int j = 0; j < s.ny; ++j )
		for( int i = 0; i < s.nx; ++i )
			b2 += s.Bx( i, j ) * s.Bx( i, j ) + s.By( i, j ) * s.By( i, j ) + s.Bz( i, j ) * s.Bz( i, j );
	const double guard = 0.1 * std::sqrt( b2 / ( static_cast< double >( s.nx ) * s.ny ) );
	maxRatio = 0.0;
	double sum = 0.0;
	int n = 0;
	for( int j = 1; j + 1 < s.ny; ++j )
		for( int i = 1; i + 1 < s.nx; ++i )
		{
			const double div = 0.5 * ( ( s.Bx( i + 1, j ) - s.Bx( i - 1, j ) ) + ( s.By( i, j + 1 ) - s.By( i, j - 1 ) ) );
			const double mag = std::sqrt( s.Bx( i, j ) * s.Bx( i, j ) + s.By( i, j ) * s.By( i, j ) + s.Bz( i, j ) * s.Bz( i, j ) );
			const double r   = std::abs( div ) / std::max( mag, guard );
			maxRatio         = std::max( maxRatio, r );
			sum += r * r;
			++n;
		}
	rmsRatio = std::sqrt( sum / std::max( n, 1 ) );
}
} // namespace

int RunDivB( const Perturb& perturb )
{
	Say( "\n=== divb: |div B| dx / |B| on a turbulent run: strong curvature, the drive, the coils turning\n" );
	for( int cells : { 128, 256 } )
	{
		TestModel model;
		model.glm = !perturb.glmOff;
		Rig rig;
		rig.plugin.SetModelForTest( model );
		rig.Init( 320, 180 );
		rig.Set( PT_DETAIL, DetailParam( cells ) );
		rig.Set( PT_CURVATURE, CurvatureParam( 1.0 ) );
		rig.Set( PT_DRIVE, 0.5f );
		rig.Set( PT_COIL_SPIN, 0.8f );
		double worst = 0.0, worstRms = 0.0;
		bool finite  = true;
		for( int chunk = 0; chunk < 6; ++chunk )
		{
			rig.Render( 50 );
			const Snapshot s = Snapshot::Take( rig.plugin );
			finite           = finite && s.Finite();
			double mx = 0.0, rms = 0.0;
			DivergenceOf( s, mx, rms );
			worst    = std::max( worst, mx );
			worstRms = std::max( worstRms, rms );
		}
		//The bounds, in cells: a field whose divergence is O(1) of |B| / dx is
		//a monopole a cell across -- the scheme has lost the field. GLM keeps
		//it at the level of the truncation error the cleaning waves carry
		//away: a tenth of |B| per cell at the very worst cell (a shock or a
		//finger a few cells wide) and a hundredth in the rms.
		Check( finite && worst < 0.1 && worstRms < 0.01,
		       fmt( "Detail %d%s: max %.4f, rms %.5f over 300 frames (bounds 0.1 and 0.01)", cells,
		            perturb.glmOff ? ", GLM OFF" : "", worst, worstRms ) );
	}
	return g_failures;
}

//===========================================================================
// --still
//===========================================================================
int RunStill( const Perturb& perturb )
{
	Say( "\n=== still: Mix 0 is the identity, bit for bit\n" );
	for( int raster : { 0, 1 } )
	{
		const int w = raster ? 1280 : 480, h = raster ? 720 : 270;
		//Every value a float can hold in 0..1, not just 8-bit steps.
		Floats card( static_cast< size_t >( w ) * h * 4 );
		uint32_t x = 0x12345678u;
		for( float& v : card )
		{
			x ^= x << 13;
			x ^= x >> 17;
			x ^= x << 5;
			v = static_cast< float >( x & 0xffffff ) / 16777216.0f;
		}
		Rig rig;
		rig.Init( w, h, &card );
		rig.Set( PT_MIX, 0.0f );
		rig.Render( 30 );
		Floats out = rig.Output();
		Floats expect = card;
		if( perturb.stillUlp )
			expect[ 4 * ( w * ( h / 2 ) + w / 2 ) ] = std::nextafter( expect[ 4 * ( w * ( h / 2 ) + w / 2 ) ], 2.0f );
		size_t differ = 0;
		for( size_t i = 0; i < out.size(); ++i )
			if( std::memcmp( &out[ i ], &expect[ i ], sizeof( float ) ) != 0 )
				++differ;
		Check( differ == 0, fmt( "%dx%d, 30 frames of plasma running underneath: %zu of %zu values differ from the input",
		                         w, h, differ, out.size() ) );
	}
	return g_failures;
}

//===========================================================================
// --state
//===========================================================================
int RunState( const Perturb& )
{
	Say( "\n=== state: what the host hands over is what it gets back\n" );
	Rig rig;
	if( !rig.Init( 320, 180 ) )
		return 1;
	rig.Set( PT_FIELD_LINES, 0.5f );
	GLuint hostArray = 0;
	glGenVertexArrays( 1, &hostArray );
	int problems = 0;
	std::string what;
	for( int frame = 0; frame < 3; ++frame )
	{
		glBindFramebuffer( GL_FRAMEBUFFER, rig.outputFBO );
		glViewport( 7, 5, 300, 170 );
		glBindVertexArray( hostArray );
		glEnable( GL_BLEND );
		glBlendFuncSeparate( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO );
		glClearColor( 0.2f, 0.3f, 0.4f, 0.5f );
		glEnable( GL_SCISSOR_TEST );
		glScissor( 0, 0, 320, 180 );
		glActiveTexture( GL_TEXTURE0 );
		glUseProgram( 0 );
		rig.plugin.SetTime( frame / 60.0 );
		if( rig.plugin.ProcessOpenGL( &rig.process ) != FF_SUCCESS )
			return 1;
		GLint viewport[ 4 ] = {}, array = 0, program = 0, unit = 0, fbo = 0, src = 0, dst = 0;
		GLfloat clear[ 4 ] = {};
		glGetIntegerv( GL_VIEWPORT, viewport );
		glGetIntegerv( GL_VERTEX_ARRAY_BINDING, &array );
		glGetIntegerv( GL_CURRENT_PROGRAM, &program );
		glGetIntegerv( GL_ACTIVE_TEXTURE, &unit );
		glGetIntegerv( GL_FRAMEBUFFER_BINDING, &fbo );
		glGetIntegerv( GL_BLEND_SRC_RGB, &src );
		glGetIntegerv( GL_BLEND_DST_RGB, &dst );
		glGetFloatv( GL_COLOR_CLEAR_VALUE, clear );
		auto expect = [ & ]( bool ok, const char* name ) {
			if( !ok )
			{
				++problems;
				what += std::string( " " ) + name;
			}
		};
		expect( viewport[ 0 ] == 7 && viewport[ 1 ] == 5 && viewport[ 2 ] == 300 && viewport[ 3 ] == 170, "viewport" );
		expect( array == static_cast< GLint >( hostArray ), "vertex-array" );
		expect( program == 0, "program" );
		expect( unit == GL_TEXTURE0, "active-unit" );
		expect( fbo == static_cast< GLint >( rig.outputFBO ), "framebuffer" );
		expect( glIsEnabled( GL_BLEND ) && src == GL_SRC_ALPHA && dst == GL_ONE_MINUS_SRC_ALPHA, "blend" );
		expect( glIsEnabled( GL_SCISSOR_TEST ), "scissor" );
		expect( clear[ 0 ] == 0.2f && clear[ 1 ] == 0.3f && clear[ 2 ] == 0.4f && clear[ 3 ] == 0.5f, "clear-colour" );
		for( int u = 0; u < 16; ++u )
		{
			GLint bound = 0;
			glActiveTexture( static_cast< GLenum >( GL_TEXTURE0 + u ) );
			glGetIntegerv( GL_TEXTURE_BINDING_2D, &bound );
			expect( bound == 0, "texture-unit" );
		}
		glActiveTexture( GL_TEXTURE0 );
	}
	glDisable( GL_SCISSOR_TEST );
	glDisable( GL_BLEND );
	glBindVertexArray( 0 );
	glDeleteVertexArrays( 1, &hostArray );
	Check( problems == 0, fmt( "three frames: viewport, vertex array, program, active unit, framebuffer, blend, scissor, "
	                           "clear colour and sixteen texture units all as the host left them (%d wrong:%s)",
	                           problems, what.empty() ? " none" : what.c_str() ) );
	return g_failures;
}

//===========================================================================
// --glow
//===========================================================================
int RunGlow( const Perturb& perturb )
{
	Say( "\n=== glow: the glare moves light out of the core and makes none\n" );
	//On the grid, where the light is made: every Gaussian stage of the glare,
	//summed over its (reduced) grid and scaled back up by the reduction, must
	//hold exactly the emission's total. A normalised kernel reflected at a
	//half-sample mirror is a symmetric operator, so this is exact but for
	//rounding: each output is a sum of up to 2 x 116 + 1 products, over two
	//passes per stage and three cascaded stages -- 6 x 233 u is the bound.
	{
		Rig rig;
		rig.Init( 640, 360 );
		rig.Set( PT_GLOW, 0.8f );
		rig.Set( PT_DETAIL, DetailParam( 512 ) );
		rig.Render( 120 );
		const Grid& g  = rig.plugin.CurrentGrid();
		const Floats e = ReadTexture( rig.plugin.EmissionTextureID(), g.nx, g.ny );
		double total   = 0.0;
		for( size_t i = 0; i < e.size(); i += 4 )
			total += e[ i ] + e[ i + 1 ] + e[ i + 2 ];
		const int gw = rig.plugin.GlowWidth(), gh = rig.plugin.GlowHeight();
		const double scale = static_cast< double >( g.nx ) * g.ny / ( static_cast< double >( gw ) * gh );
		const double bound = 6.0 * 233.0 * std::ldexp( 1.0, -24 );
		for( int stage = 0; stage < kGlowStages; ++stage )
		{
			const Floats glow = ReadTexture( rig.plugin.GlowTextureID( stage ), gw, gh );
			double sum        = 0.0;
			for( size_t i = 0; i < glow.size(); i += 4 )
				sum += glow[ i ] + glow[ i + 1 ] + glow[ i + 2 ];
			sum *= scale;
			Check( std::abs( sum / total - 1.0 ) < bound,
			       fmt( "stage %d (sigma %.3f of the frame height, on a %dx%d copy of the %dx%d grid): %.8f of the emission "
			            "(bound 1 +- %.1e)",
			            stage, kGlowSigma[ stage ], gw, gh, g.nx, g.ny, sum / total, bound ) );
		}
	}
	//On the raster: the whole picture's light with Glow 0.6 and Glow 0,
	//exposure low enough that nothing reaches the shoulder, at a raster twice
	//the grid (where bilinear resampling keeps a sum exactly). A bloom that
	//added its glare would read 1.6.
	{
		double totals[ 2 ] = {};
		for( int k = 0; k < 2; ++k )
		{
			TestModel model;
			model.additiveGlow = perturb.additiveGlow;
			Rig rig;
			rig.plugin.SetModelForTest( model );
			rig.Init( 512, 512 );
			rig.Set( PT_DETAIL, DetailParam( 256 ) );
			rig.Set( PT_EXPOSURE, 0.0f );
			rig.Set( PT_SPEED, 0.0f );
			rig.Set( PT_GLOW, k == 0 ? 0.0f : 0.667f );
			rig.Render( 2 );
			const Floats out = rig.Output();
			double peak      = 0.0;
			for( size_t i = 0; i < out.size(); i += 4 )
			{
				totals[ k ] += out[ i ] + out[ i + 1 ] + out[ i + 2 ];
				peak = std::max( { peak, static_cast< double >( out[ i ] ), static_cast< double >( out[ i + 1 ] ),
				                   static_cast< double >( out[ i + 2 ] ) } );
			}
			if( peak > 0.8 )
				Check( false, fmt( "the raster check needs every pixel below the shoulder's knee; peak %.3f", peak ) );
		}
		Check( std::abs( totals[ 1 ] / totals[ 0 ] - 1.0 ) < 1e-4,
		       fmt( "the raster's total light with Glow 0.6 is %.6f of it with none (bound 1 +- 1e-4)",
		            totals[ 1 ] / totals[ 0 ] ) );
	}
	return g_failures;
}

#define STUB( name )                                  \
	int name( const Perturb& )                        \
	{                                                 \
		Check( false, #name " not written yet" );     \
		return 1;                                     \
	}
STUB( RunBalance ) STUB( RunRT ) STUB( RunCusp ) STUB( RunFrozen )
STUB( RunQuench ) STUB( RunResist ) STUB( RunFloors ) 
STUB( RunMutation )

//===========================================================================
// --negative
//===========================================================================
int RunNegative()
{
	struct Case
	{
		const char* name;
		int ( *check )( const Perturb& );
		Perturb perturb;
		const char* what;
	};
	std::vector< Case > cases;
	auto add = [ & ]( const char* name, int ( *check )( const Perturb& ), const char* what,
	                  const std::function< void( Perturb& ) >& set ) {
		Perturb p;
		set( p );
		cases.push_back( { name, check, p, what } );
	};
	add( "briowu", RunBrioWu, "run HLL where HLLD is claimed", []( Perturb& p ) { p.hll = true; } );
	add( "alfven", RunAlfven, "expect the Alfven speed B / rho", []( Perturb& p ) { p.alfvenOverRho = true; } );
	add( "conserve", RunConserve, "run with cooling, a real energy sink", []( Perturb& p ) { p.conserveCooling = true; } );
	add( "divb", RunDivB, "run with GLM off", []( Perturb& p ) { p.glmOff = true; } );
	add( "still", RunStill, "expect one ulp more in one pixel", []( Perturb& p ) { p.stillUlp = true; } );
	add( "glow", RunGlow, "run the glare as an additive bloom", []( Perturb& p ) { p.additiveGlow = true; } );
	if( const char* only = std::getenv( "CT_NEGATIVE" ) )
		cases.erase( std::remove_if( cases.begin(), cases.end(), [ & ]( const Case& c ) { return std::string( c.name ) != only; } ),
		             cases.end() );

	int unfalsifiable = 0;
	for( const Case& c : cases )
	{
		Say( "\n=== negative control: %s -- %s\n", c.name, c.what );
		const int before = g_failures;
		g_failures       = 0;
		c.check( c.perturb );
		const int observed = g_failures;
		g_failures         = before;
		if( observed > 0 )
			Say( "  ok    %s failed %d check%s, as it must\n", c.name, observed, observed == 1 ? "" : "s" );
		else
		{
			Say( "  FAIL  %s PASSED against a wrong model -- it cannot fail, so it is not a check\n", c.name );
			++unfalsifiable;
		}
	}
	Say( "\nnegative controls: %zu wrong models, %d of them undetected\n", cases.size(), unfalsifiable );
	Say( "\n  %s\n", unfalsifiable == 0 ? "PASS" : "FAIL" );
	return unfalsifiable == 0 ? 0 : 1;
}


//===========================================================================
// --bench
//===========================================================================
int RunBench()
{
	Say( "\n=== bench: ms per frame, the default look, 60 frames after 30 of warm-up\n" );
	struct Size
	{
		int w, h;
		const char* name;
	};
	const Size sizes[] = { { 1280, 720, "720p" }, { 1920, 1080, "1080p" }, { 3840, 2160, "4K" } };
	for( int cells : { 128, 256, 512, 1024 } )
		for( const Size& size : sizes )
		{
			Rig rig;
			if( !rig.Init( size.w, size.h ) )
				return 1;
			rig.Set( PT_DETAIL, DetailParam( cells ) );
			if( !rig.Render( 30 ) )
				return 1;
			glFinish();
			const long long before = rig.plugin.SubstepsTaken();
			const auto start       = std::chrono::steady_clock::now();
			if( !rig.Render( 60 ) )
				return 1;
			glFinish();
			const auto end  = std::chrono::steady_clock::now();
			const double ms = std::chrono::duration< double, std::milli >( end - start ).count() / 60.0;
			const Grid& g   = rig.plugin.CurrentGrid();
			Say( "  Detail %-4d %-6s %7.2f ms/frame  (%5.1f%% of 60 fps; grid %dx%d, %.1f substeps/frame, %d frames capped)\n",
			     cells, size.name, ms, 100.0 * ms / ( 1000.0 / 60.0 ), g.nx, g.ny,
			     ( rig.plugin.SubstepsTaken() - before ) / 60.0, rig.plugin.CappedFrames() );
		}
	return 0;
}


} // namespace cttest
namespace cttest
{
//A development aid: the coils' vacuum field in a uniform background should be
//an equilibrium. Prints how far from one it drifts.
int RunEquilibrium( const Perturb& )
{
	for( int boundary : { 0, 1 } )
		for( int glm : { 1, 0 } )
			for( int poles : { 0, 4, 6 } )
			{
				TestModel model;
				model.boundary = boundary;
				model.glm      = glm == 1;
				Rig rig;
				rig.plugin.SetModelForTest( model );
				rig.Init( 256, 144 );
				rig.Quiet();
				rig.Set( PT_POLES, PolesParam( poles ) );
				rig.Set( PT_SPEED, 0.0f );
				rig.Render( 1 );
				const Grid& g = rig.plugin.CurrentGrid();
				const Coils& coils = rig.plugin.CurrentCoils();
				StateBuilder sb( g.nx, g.ny, kGamma );
				for( int j = 0; j < g.ny; ++j )
					for( int i = 0; i < g.nx; ++i )
					{
						double bx, by, bz;
						VacuumField( coils, ( i + 0.5 ) * g.dx, ( j + 0.5 ) * g.dx, bx, by, bz );
						sb.Set( i, j, kBackgroundDensity, 0, 0, 0, kBackgroundPressure, bx, by, bz );
					}
				sb.Load( rig.plugin );
				const int steps = rig.plugin.StepForTest( std::getenv( "CT_T" ) ? std::atof( std::getenv( "CT_T" ) ) : 0.05 );
				const Snapshot s = Snapshot::Take( rig.plugin );
				double vmax = 0, pmin = 1e30, pmax = 0;
				int wi = 0, wj = 0;
				for( int j = 0; j < s.ny; ++j )
					for( int i = 0; i < s.nx; ++i )
					{
						const double v = std::hypot( s.U( i, j ), s.V( i, j ) );
						if( v > vmax )
						{
							vmax = v;
							wi   = i;
							wj   = j;
						}
						pmin = std::min( pmin, s.P( i, j ) );
						pmax = std::max( pmax, s.P( i, j ) );
					}
				Say( "  boundary %d glm %d poles %d: %d steps, |v|max %.3g at (%d,%d) of %dx%d, p [%.3g %.3g], floors %.0f\n",
				     boundary, glm, poles, steps, vmax, wi, wj, s.nx, s.ny, pmin, pmax, rig.plugin.FloorHits() );
			}
	return 0;
}
}
