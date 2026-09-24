#include "Checks.h"

#include "Presets.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <chrono>
#include <cstdio>
#include <functional>
#include <limits>
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
		int w = raster ? 1280 : 480, h = raster ? 720 : 270;
		ChooseRaster( w, h );
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
		const int ew = rig.plugin.EmissionWidth(), eh = rig.plugin.EmissionHeight();
		const Floats e = ReadTexture( rig.plugin.EmissionTextureID(), ew, eh );
		double total   = 0.0;
		for( size_t i = 0; i < e.size(); i += 4 )
			total += e[ i ] + e[ i + 1 ] + e[ i + 2 ];
		const int gw = rig.plugin.GlowWidth(), gh = rig.plugin.GlowHeight();
		const double scale = static_cast< double >( ew ) * eh / ( static_cast< double >( gw ) * gh );
		const double bound = 6.0 * 233.0 * std::ldexp( 1.0, -24 );
		for( int stage = 0; stage < kGlowStages; ++stage )
		{
			const Floats glow = ReadTexture( rig.plugin.GlowTextureID( stage ), gw, gh );
			double sum        = 0.0;
			for( size_t i = 0; i < glow.size(); i += 4 )
				sum += glow[ i ] + glow[ i + 1 ] + glow[ i + 2 ];
			sum *= scale;
			Check( std::abs( sum / total - 1.0 ) < bound,
			       fmt( "stage %d (sigma %.3f of the frame height, on a %dx%d copy of the %dx%d frame): %.8f of the emission "
			            "(bound 1 +- %.1e)",
			            stage, kGlowSigma[ stage ], gw, gh, ew, eh, sum / total, bound ) );
		}
	}
	//On the raster: the whole picture's light with Glow 0.6 and Glow 0,
	//exposure low enough that nothing reaches the shoulder. A bloom that
	//added its glare would read 1.6.
	//
	//The raster samples each grid bilinearly, and a resample conserves a sum
	//only where every cell collects the same weight from the pixels. At a
	//raster exactly twice the grid it does (but at the clamped border), and
	//the check was written there; at 320x180 it does not. So the bound is
	//the resampler's own, computed on the CPU for these sizes and weighted
	//by where the light actually is: for each texture the pixels' summed
	//weight on each cell, W, against its mean, plus 2^-9 a pixel wherever the
	//filter's fraction is not a multiple of 1/256 -- GL leaves the filter's
	//sub-texel precision to the implementation, and 8 bits is the common
	//floor. Ten float roundings a pixel, correlated at worst, on each total.
	{
		struct Axis
		{
			std::vector< double > w, q;
		};
		auto axis = []( int g, int r ) {
			Axis out { std::vector< double >( g, 0.0 ), std::vector< double >( g, 0.0 ) };
			for( int p = 0; p < r; ++p )
			{
				const double s = ( p + 0.5 ) * g / r - 0.5, f = s - std::floor( s );
				const int i0 = static_cast< int >( std::floor( s ) );
				const int c0 = std::clamp( i0, 0, g - 1 ), c1 = std::clamp( i0 + 1, 0, g - 1 );
				const double q = f * 256.0 == std::floor( f * 256.0 ) ? 0.0 : std::ldexp( 1.0, -9 );
				out.w[ c0 ] += 1.0 - f;
				out.w[ c1 ] += f;
				out.q[ c0 ] += q;
				out.q[ c1 ] += q;
			}
			return out;
		};
		//The relative error of one texture's resampled total, for its light.
		auto ripple = [ & ]( const Floats& light, int gw, int gh, int rw, int rh ) {
			const Axis x = axis( gw, rw ), y = axis( gh, rh );
			const double mean = static_cast< double >( rw ) * rh / ( static_cast< double >( gw ) * gh );
			double err = 0.0, total = 0.0;
			for( int j = 0; j < gh; ++j )
				for( int i = 0; i < gw; ++i )
				{
					const size_t o = 4 * ( static_cast< size_t >( j ) * gw + i );
					const double l = light[ o ] + light[ o + 1 ] + light[ o + 2 ];
					const double w = x.w[ i ] * y.w[ j ];
					const double q = x.w[ i ] * y.q[ j ] + x.q[ i ] * y.w[ j ] + x.q[ i ] * y.q[ j ];
					err += ( std::abs( w - mean ) + q ) * l;
					total += l;
				}
			return total > 0.0 ? err / ( mean * total ) : 0.0;
		};
		double totals[ 2 ] = {}, spread[ 2 ] = {};
		int rw = 512, rh = 512;
		ChooseRaster( rw, rh );
		//The grid no bigger than the raster: sampled DOWN, a bilinear filter
		//skips cells outright and no sum survives at all (the ripple is ~1).
		const int cells = std::min( rw, rh ) >= 512 ? 256 : 128;
		for( int k = 0; k < 2; ++k )
		{
			TestModel model;
			model.additiveGlow = perturb.additiveGlow;
			Rig rig;
			rig.plugin.SetModelForTest( model );
			rig.Init( 512, 512 );
			rig.Set( PT_DETAIL, DetailParam( cells ) );
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
			const int ew = rig.plugin.EmissionWidth(), eh = rig.plugin.EmissionHeight();
			spread[ k ] = ripple( ReadTexture( rig.plugin.EmissionTextureID(), ew, eh ), ew, eh, rw, rh );
			if( k == 1 )
				for( int stage = 0; stage < kGlowStages; ++stage )
				{
					const int gw = rig.plugin.GlowWidth(), gh = rig.plugin.GlowHeight();
					spread[ k ] += ripple( ReadTexture( rig.plugin.GlowTextureID( stage ), gw, gh ), gw, gh, rw, rh );
				}
		}
		const double rounding = 2.0 * 10.0 * std::ldexp( 1.0, -24 );
		const double bound    = ( spread[ 1 ] + spread[ 0 ] ) / ( 1.0 - spread[ 0 ] ) + rounding;
		Check( std::abs( totals[ 1 ] / totals[ 0 ] - 1.0 ) < bound,
		       fmt( "at %dx%d, Detail %d, the raster's total light with Glow 0.6 is %.6f of it with none (bound 1 +- %.1e: "
		            "the resampler's own ripple where the light is, %.1e and %.1e, and rounding)",
		            rw, rh, cells, totals[ 1 ] / totals[ 0 ], bound, spread[ 1 ], spread[ 0 ] ) );
	}
	return g_failures;
}


//===========================================================================
// Radial profiles about a centre, azimuthally averaged, in rings of one cell.
//===========================================================================
namespace
{
struct Profile
{
	std::vector< double > r, p, bz, chi, beta, rho, rhoV2;
};

Profile Radial( const Snapshot& s, double cx, double cy, double rMax )
{
	const int bins = static_cast< int >( rMax / s.dx );
	Profile out;
	std::vector< double > n( bins, 0.0 ), p( bins, 0.0 ), bz( bins, 0.0 ), chi( bins, 0.0 ), beta( bins, 0.0 ),
		rho( bins, 0.0 ), kin( bins, 0.0 );
	for( int j = 0; j < s.ny; ++j )
		for( int i = 0; i < s.nx; ++i )
		{
			const double x = ( i + 0.5 ) * s.dx - cx, y = ( j + 0.5 ) * s.dx - cy;
			const int b    = static_cast< int >( std::sqrt( x * x + y * y ) / s.dx );
			if( b >= bins )
				continue;
			const double pp = s.P( i, j );
			const double b2 = s.Bx( i, j ) * s.Bx( i, j ) + s.By( i, j ) * s.By( i, j ) + s.Bz( i, j ) * s.Bz( i, j );
			n[ b ] += 1.0;
			p[ b ] += pp;
			bz[ b ] += s.Bz( i, j );
			chi[ b ] += s.Tracer( i, j, 3 );
			beta[ b ] += 2.0 * pp / std::max( b2, 1e-30 );
			rho[ b ] += s.Rho( i, j );
			kin[ b ] += s.Rho( i, j ) * ( s.U( i, j ) * s.U( i, j ) + s.V( i, j ) * s.V( i, j ) );
		}
	for( int b = 0; b < bins; ++b )
	{
		if( n[ b ] == 0.0 )
			continue;
		out.r.push_back( ( b + 0.5 ) * s.dx );
		out.p.push_back( p[ b ] / n[ b ] );
		out.bz.push_back( bz[ b ] / n[ b ] );
		out.chi.push_back( chi[ b ] / n[ b ] );
		out.beta.push_back( beta[ b ] / n[ b ] );
		out.rho.push_back( rho[ b ] / n[ b ] );
		out.rhoV2.push_back( kin[ b ] / n[ b ] );
	}
	return out;
}

double MeanOver( const std::vector< double >& r, const std::vector< double >& v, double lo, double hi )
{
	double sum = 0.0, n = 0.0;
	for( size_t k = 0; k < r.size(); ++k )
		if( r[ k ] >= lo && r[ k ] <= hi )
		{
			sum += v[ k ];
			n += 1.0;
		}
	return n > 0.0 ? sum / n : 0.0;
}

/// Where v first falls through `level` going outwards, interpolated.
double Crossing( const std::vector< double >& r, const std::vector< double >& v, double level )
{
	for( size_t k = 0; k + 1 < r.size(); ++k )
		if( ( v[ k ] - level ) * ( v[ k + 1 ] - level ) <= 0.0 && v[ k ] != v[ k + 1 ] )
			return r[ k ] + ( r[ k + 1 ] - r[ k ] ) * ( v[ k ] - level ) / ( v[ k ] - v[ k + 1 ] );
	return -1.0;
}
} // namespace

//===========================================================================
// --balance
//===========================================================================
int RunBalance( const Perturb& perturb )
{
	Say( "\n=== balance: the diamagnetic bubble -- guide field only, no curvature, ideal\n" );
	for( int cells : { 256, 512 } )
	{
		Rig rig;
		rig.Init( 256, 256 );
		rig.Quiet();
		rig.Set( PT_DETAIL, DetailParam( cells ) );
		rig.Set( PT_SPEED, 0.0f );
		rig.Set( PT_BOUNDARY, 0.0f );//Open: the fast waves the ringing sends out leave
		rig.Set( PT_POLES, PolesParam( 0 ) );
		rig.Set( PT_FIELD, FieldParam( 1.0 ) );
		rig.Set( PT_GUIDE_FIELD, GuideParam( 1.0 ) );
		rig.Set( PT_PROFILE, 1.0f );
		rig.Set( PT_BALL_SIZE, BallSizeParam( 0.15 ) );
		rig.Set( PT_TEMPERATURE, TemperatureParam( 1.2 ) );
		if( perturb.balanceEta > 0.0 )
			rig.Set( PT_RESISTIVITY, ResistivityParam( perturb.balanceEta ) );
		rig.Render( 1 );
		const Grid& g   = rig.plugin.CurrentGrid();
		const double cx = 0.5 * g.lx, cy = 0.5 * g.ly;

		//The ball's mass, rho chi summed: the scheme conserves it exactly (its
		//fluxes telescope), and nothing reaches the Open margin in this run.
		auto ballMass = []( const Snapshot& s ) {
			double m = 0.0;
			for( int j = 0; j < s.ny; ++j )
				for( int i = 0; i < s.nx; ++i )
					m += s.Rho( i, j ) * s.Tracer( i, j, 3 ) * s.dx * s.dx;
			return m;
		};
		//The core as laid down: a top hat, so everything inside half its
		//radius is the ball's own rho0, p0, Bz0.
		const Snapshot s0   = Snapshot::Take( rig.plugin );
		const double mass0  = ballMass( s0 );
		const Profile pr0   = Radial( s0, cx, cy, 0.48 );
		const double rho0   = MeanOver( pr0.r, pr0.rho, 0.0, 0.075 );
		const double p0     = MeanOver( pr0.r, pr0.p, 0.0, 0.075 );
		const double bz0    = MeanOver( pr0.r, pr0.bz, 0.0, 0.075 );
		const double r0     = std::sqrt( mass0 / ( kPi * rho0 ) );

		const int steps  = rig.plugin.StepForTest( 8.0 );
		const Snapshot s = Snapshot::Take( rig.plugin );

		//The edge is the radius of the ball's VOLUME: its mass over its core
		//density. At pressure balance the ball's plasma, in the core and in
		//any finger alike, sits on the core's adiabat with the core's Bz / rho,
		//so it has the core's density everywhere, and M / rho_core is its
		//volume whatever shape the edge has taken. (The marker's half-contour
		//is not: the ringing's first deceleration is Rayleigh-Taylor unstable
		//-- a dense ball, a tenuous background, k across B so nothing holds it
		//-- and in the mixed layer a cell that is half ball by MASS is only
		//rho_out / ( rho_in + rho_out ) ~ 13% ball by volume. At Detail 512 the
		//layer is resolved well enough to finger, and that contour sat 4.5
		//cells out. It is printed below, not asserted.)
		const double mass = ballMass( s );
		const Profile pr  = Radial( s, cx, cy, 0.48 );
		double halfArea   = 0.0;
		for( int j = 0; j < s.ny; ++j )
			for( int i = 0; i < s.nx; ++i )
				if( s.Tracer( i, j, 3 ) > 0.5 )
					halfArea += s.dx * s.dx;
		const double rhoCore = MeanOver( pr.r, pr.rho, 0.0, 0.5 * r0 );
		const double edge    = std::sqrt( mass / ( kPi * rhoCore ) );
		const double pIn     = MeanOver( pr.r, pr.p, 0.0, 0.5 * edge );
		const double bIn     = MeanOver( pr.r, pr.bz, 0.0, 0.5 * edge );
		const double pOut    = MeanOver( pr.r, pr.p, 1.5 * edge, 2.0 * edge );
		const double bOut    = MeanOver( pr.r, pr.bz, 1.5 * edge, 2.0 * edge );
		const double totalIn = pIn + 0.5 * bIn * bIn, totalOut = pOut + 0.5 * bOut * bOut;
		//What ringing is left: the momentum equation balances any total-
		//pressure difference against rho dv/dt ~ rho v^2 / L, so the largest
		//rho v^2 anywhere bounds it. That plus half a percent for the edge's
		//own numerical width is the tolerance, relative to the total.
		double dyn = 0.0;
		for( double v : pr.rhoV2 )
			dyn = std::max( dyn, v );
		const double tol = dyn / totalOut + 0.005;
		Check( std::abs( totalIn / totalOut - 1.0 ) < tol,
		       fmt( "Detail %d, t = 8 (%d substeps): p + B^2/2 inside %.5f, outside %.5f (%.3f%%; bound %.3f%%)", cells,
		            steps, totalIn, totalOut, 100 * std::abs( totalIn / totalOut - 1.0 ), 100 * tol ) );
		const double predicted = perturb.balanceNoTwo ? std::sqrt( std::max( bOut * bOut + 2.0 * pOut - pIn, 0.0 ) )
		                                              : std::sqrt( std::max( bOut * bOut + 2.0 * pOut - 2.0 * pIn, 0.0 ) );
		Check( std::abs( bIn / predicted - 1.0 ) < tol,
		       fmt( "Detail %d: Bz inside %.5f against sqrt( B0^2 + 2 p_out - 2 p_in ) = %.5f%s (%.3f%%; bound %.3f%%)",
		            cells, bIn, predicted, perturb.balanceNoTwo ? " (without the 2)" : "",
		            100 * std::abs( bIn / predicted - 1.0 ), 100 * tol ) );

		//Where flux conservation puts the edge. In 2.5-D Bz / rho is carried
		//with the plasma, so the ball, having expanded by x in area, holds
		//Bz0 / x: x = Bz0 / Bz_core, and the edge is r0 sqrt( x ), r0 being the
		//ball's volume radius at ignition. (Not an adiabatic prediction: the
		//ring-down's compressions converge on the axis and heat the core --
		//its entropy is printed below -- so p alone would not say where the
		//edge goes. The flux does, and that is the claim.)
		//
		//The tolerance is one cell. A fast wave across B moves Bz and rho in
		//proportion, so the ringing leaves Bz / rho untouched, and the only
		//plasma the core does not describe is what the ignition's tanh edge
		//(half-width 0.75 cells) laid down part-mixed with the background: less
		//than a cell's width of the ball, so less than a cell of radius.
		const double x          = bz0 / bIn;
		const double rPredicted = r0 * std::sqrt( x );
		Check( std::abs( edge - rPredicted ) < s.dx,
		       fmt( "Detail %d: the edge (the ball's volume, M / rho_core) at r = %.5f; flux conservation, Bz0 / Bz_core = "
		            "%.4f, puts it at %.5f%s (%.2f cells; bound 1)",
		            cells, edge, x, rPredicted, perturb.balanceEta > 0.0 ? " (WITH RESISTIVITY: the field diffuses in)" : "",
		            std::abs( edge - rPredicted ) / s.dx ) );
		Say( "  note  Detail %d: the marker's half-contour encloses r = %.4f, %.1f cells outside the edge -- the mixed "
		     "layer, which that contour counts at ~13%% ball by volume\n",
		     cells, std::sqrt( halfArea / kPi ), ( std::sqrt( halfArea / kPi ) - edge ) / s.dx );
		//The beta = 1 contour lies in the edge layer: between where the ball
		//marker has fallen to 0.9 and to 0.1 (plus a cell each side), with
		//beta above 1 in the core and below 1 outside. (It is not at the
		//marker's half: with beta_in ~ 1.4, 2p = B^2 where p has fallen only
		//a tenth of the way across the numerically widened edge.)
		const double beta1  = Crossing( pr.r, pr.beta, 1.0 );
		const double inner = Crossing( pr.r, pr.chi, 0.9 ), outer = Crossing( pr.r, pr.chi, 0.1 );
		const double betaIn = 2.0 * pIn / ( bIn * bIn ), betaOut = 2.0 * pOut / ( bOut * bOut );
		Check( beta1 > inner - s.dx && beta1 < outer + s.dx && betaIn > 1.0 && betaOut < 1.0,
		       fmt( "Detail %d: beta = 1 at r = %.4f, inside the edge layer %.4f..%.4f; beta %.2f inside, %.3f outside",
		            cells, beta1, inner, outer, betaIn, betaOut ) );
		const double kIn = pIn / std::pow( rhoCore, kGamma ), kStart = p0 / std::pow( rho0, kGamma );
		Say( "  note  Detail %d: the core's entropy p / rho^gamma %.5f, %.5f at ignition (%+.3f%%)\n", cells, kIn, kStart,
		     100.0 * ( kIn / kStart - 1.0 ) );
	}
	return g_failures;
}

//===========================================================================
// --rt
//===========================================================================
int RunRT( const Perturb& perturb )
{
	Say( "\n=== rt: magnetic Rayleigh-Taylor, a planar slab, one seeded mode\n" );
	const double g = 1.0, rhoH = 2.0, rhoL = 1.0, A = ( rhoH - rhoL ) / ( rhoH + rhoL ), p0 = 20.0;
	const double k = 2.0 * kPi * 2.0;//two wavelengths across the unit box
	const double hydro = g * k * A;
	struct Case
	{
		const char* name;
		double bx, bz;
		int cells;
	};
	const double bHalf = std::sqrt( hydro * ( rhoH + rhoL ) / ( 4.0 * k * k ) ); //gamma^2 = gkA / 2
	const double bStop = std::sqrt( hydro * ( rhoH + rhoL ) / ( k * k ) );       //gamma^2 = -gkA
	const Case cases[] = { { "B perpendicular to k (guide field)", 0.0, 1.0, 256 },
		                   { "B perpendicular to k (guide field)", 0.0, 1.0, 512 },
		                   { "B along the interface, k < k_c", bHalf, 0.0, 256 },
		                   { "B along the interface, k < k_c", bHalf, 0.0, 512 },
		                   { "B along the interface, k > k_c", bStop, 0.0, 256 } };
	for( const Case& c : cases )
	{
		TestModel model;
		model.boundary       = 4;
		model.gravityCentre  = true;
		model.gravityX       = 0.5;
		model.gravityY       = -1.0e4;//g along +y, uniform to 1e-9 over the box
		model.gravityCore    = 0.0;
		model.uniformGravity = true;
		//Two slabs, identical but for the seed. The discrete hydrostatic
		//balance is not exact across the interface, so an unseeded slab
		//makes small flows of its own at every wavelength; in the linear
		//regime the seeded slab is those plus the mode, so the mode is the
		//DIFFERENCE of the two -- the spurious flows cancel exactly.
		Rig rigs[ 2 ];
		const double delta = 1e-3;
		Grid gr;
		for( int r = 0; r < 2; ++r )
		{
			Rig& rig = rigs[ r ];
			rig.plugin.SetModelForTest( model );
			rig.Init( 128, 128 );
			rig.Quiet();
			rig.Set( PT_DETAIL, DetailParam( c.cells ) );
			rig.Set( PT_SPEED, 0.0f );
			rig.Set( PT_POLES, PolesParam( 0 ) );
			rig.Set( PT_GUIDE_FIELD, 0.0f );
			rig.Set( PT_CURVATURE, CurvatureParam( g ) );
			rig.Render( 1 );
			gr = rig.plugin.CurrentGrid();
			StateBuilder sb( gr.nx, gr.ny, kGamma );
			const double w = gr.dx;//the interface: a tanh one cell wide
			for( int j = 0; j < gr.ny; ++j )
				for( int i = 0; i < gr.nx; ++i )
				{
					const double x = ( i + 0.5 ) * gr.dx, y = ( j + 0.5 ) * gr.dx - 0.5;
					//Heavy below (g points up, +y), in hydrostatic balance:
					//dp/dy = rho g, integrated in closed form.
					const double rho = 0.5 * ( rhoH + rhoL ) + 0.5 * ( rhoH - rhoL ) * std::tanh( -y / w );
					const double p   = p0 + g * ( 0.5 * ( rhoH + rhoL ) * y
					                            - 0.5 * ( rhoH - rhoL ) * w * std::log( std::cosh( y / w ) ) );
					//The incompressible eigenfunction, divergence-free.
					const double e = r == 0 ? std::exp( -k * std::abs( y ) ) : 0.0;
					const double v = delta * std::cos( k * x ) * e;
					const double u = delta * ( y >= 0.0 ? 1.0 : -1.0 ) * std::sin( k * x ) * e;
					sb.Set( i, j, rho, u, v, 0.0, p, c.bx, 0.0, c.bz );
				}
			sb.Load( rig.plugin );
		}
		auto amplitude = [ & ]() {
			const Snapshot s = Snapshot::Take( rigs[ 0 ].plugin );
			const Snapshot z = Snapshot::Take( rigs[ 1 ].plugin );
			double re = 0.0, im = 0.0;
			for( int j = 0; j < s.ny; ++j )
			{
				const double y = ( j + 0.5 ) * s.dx - 0.5;
				if( std::abs( y ) > 0.15 )
					continue;
				for( int i = 0; i < s.nx; ++i )
				{
					const double x = ( i + 0.5 ) * s.dx;
					const double v = s.V( i, j ) - z.V( i, j );
					re += v * std::cos( k * x );
					im += v * std::sin( k * x );
				}
			}
			return std::sqrt( re * re + im * im );
		};
		std::vector< double > t, a;
		const double a0 = amplitude();
		double peak     = a0;
		//Run for three e-folds of the theory's rate (or 2 tau_A if stable),
		//sampled 25 times.
		const double gRate = std::sqrt( std::max( perturb.curvatureSign * g * k * A - 2.0 * k * k * c.bx * c.bx / ( rhoH + rhoL ), 0.0 ) );
		const double span  = gRate > 0.0 ? 3.0 / gRate : 2.0;
		for( int q = 1; q <= 25; ++q )
		{
			rigs[ 0 ].plugin.StepForTest( span / 25.0 );
			rigs[ 1 ].plugin.StepForTest( span / 25.0 );
			t.push_back( span * q / 25.0 );
			a.push_back( amplitude() );
			peak = std::max( peak, a.back() );
		}
		//Theory, with the model's sign of g (the negative control flips it).
		const double gs = perturb.curvatureSign * g;
		const double gamma2 = gs * k * A - 2.0 * k * k * c.bx * c.bx / ( rhoH + rhoL );
		if( gamma2 > 0.0 )
		{
			//The seed is velocity with no displacement: equal parts of the
			//growing and the decaying mode, so v = v0 cosh( gamma t ), not
			//v0 e^( gamma t ). (Fitting the slope of ln v -- the first
			//version -- read gamma tanh( gamma t ): 6% low at two e-folds.)
			//gamma from each sample in the last two e-folds, acosh( v / v0 ) / t,
			//and the median of them.
			std::vector< double > rates;
			for( size_t q = 0; q < t.size(); ++q )
				if( t[ q ] >= span / 3.0 - 1e-9 && a[ q ] > a0 )
					rates.push_back( std::acosh( a[ q ] / a0 ) / t[ q ] );
			std::sort( rates.begin(), rates.end() );
			const double rate = rates.empty() ? 0.0 : rates[ rates.size() / 2 ];
			const double theory = std::sqrt( gamma2 );
			//The bound: a diffuse interface of thickness L weakens the
			//buoyancy term to gkA / ( 1 + kL ), and the scheme holds this
			//interface about two cells thick. The tension term does not
			//depend on L, so the relative error in the rate is
			//( kL / 2 ) gkA / gamma^2: twice as large where tension has
			//taken half of gkA. Plus compressibility, g / ( k c_s^2 ), and 2%
			//for the fit.
			const double cs2 = kGamma * p0 / rhoL;
			const double tol = 0.5 * k * 2.0 * gr.dx * hydro / gamma2 + g / ( k * cs2 ) + 0.02;
			Check( std::abs( rate / theory - 1.0 ) < tol,
			       fmt( "%s, Detail %d: growth %.4f against sqrt( gkA - 2 (k.B)^2 / (rho1 + rho2) ) = %.4f (%.2f%%; bound %.2f%%)",
			            c.name, c.cells, rate, theory, 100 * std::abs( rate / theory - 1.0 ), 100 * tol ) );
		}
		else
		{
			//Stable: it oscillates at sqrt( -gamma^2 ) and must not grow. The
			//unstable cases grow by e^3 over the same time.
			Check( peak < 2.0 * a0, fmt( "%s, Detail %d: gamma^2 = %.3f < 0 predicted; amplitude peaked at %.2f of its start "
			                             "(bound 2; the unstable cases reach > 20)",
			                             c.name, c.cells, gamma2, peak / a0 ) );
		}
	}
	return g_failures;
}

//===========================================================================
// --cusp
//===========================================================================
namespace
{
/// The cusp angles on a circle, from the coils alone: where the vacuum field
/// is radial (B_phi = 0), which is where its field lines lead out. Computed
/// here in double from the line currents, not from the plugin's shader.
std::vector< double > CuspAngles( const Coils& coils, double cx, double cy, double radius )
{
	std::vector< double > out;
	const int n = 7200;
	auto bphi = [ & ]( double a ) {
		double bx, by, bz;
		VacuumField( coils, cx + radius * std::cos( a ), cy + radius * std::sin( a ), bx, by, bz );
		return -bx * std::sin( a ) + by * std::cos( a );
	};
	for( int i = 0; i < n; ++i )
	{
		const double a0 = 2.0 * kPi * i / n, a1 = 2.0 * kPi * ( i + 1 ) / n;
		const double f0 = bphi( a0 ), f1 = bphi( a1 );
		if( f0 * f1 < 0.0 )
		{
			//Keep the zeros where the field points OUT or IN along the radius
			//and is strong -- the cusps; B_phi also vanishes, weakly, nowhere
			//else on a pure multipole, but check |B_r| anyway.
			const double a = a0 + ( a1 - a0 ) * f0 / ( f0 - f1 );
			double bx, by, bz;
			VacuumField( coils, cx + radius * std::cos( a ), cy + radius * std::sin( a ), bx, by, bz );
			if( std::hypot( bx, by ) > 1e-6 )
				out.push_back( a );
		}
	}
	return out;
}

double AngleDiff( double a, double b )
{
	double d = std::fmod( a - b, 2.0 * kPi );
	if( d > kPi )
		d -= 2.0 * kPi;
	if( d < -kPi )
		d += 2.0 * kPi;
	return d;
}
} // namespace

int RunCusp( const Perturb& perturb )
{
	Say( "\n=== cusp: a hot ball in an N-pole cusp leaks where the vacuum field says\n" );
	const double rc = 0.33;
	struct Case
	{
		int poles;
		double spin;
	};
	const Case cases[] = { { 4, 0.0 }, { 6, 0.0 }, { 8, 0.0 }, { 6, 0.5 } };
	for( const Case& c : cases )
	{
		Rig rig;
		rig.Init( 256, 256 );
		rig.Quiet();
		rig.Set( PT_DETAIL, DetailParam( 256 ) );
		//Open, so what leaks leaves; with the coils turning, Wall, where the
		//coils' turn is carried through the whole field at once. (With Open
		//it arrives only through the edge and the interior lags it.)
		rig.Set( PT_BOUNDARY, c.spin > 0.0 ? 1.0f : 0.0f );
		rig.Set( PT_POLES, PolesParam( c.poles ) );
		rig.Set( PT_FIELD, FieldParam( 1.0 ) );
		rig.Set( PT_GUIDE_FIELD, 0.0f );
		rig.Set( PT_BALL_SIZE, BallSizeParam( 0.12 ) );
		rig.Set( PT_TEMPERATURE, TemperatureParam( 2.0 ) );
		rig.Set( PT_COIL_SPIN, static_cast< float >( 0.5 + 0.5 * c.spin ) );
		rig.Set( PT_SPEED, ParamFromSpeed( 0.6f ) );//0.01 tau_A a frame
		const Grid& g = rig.plugin.CurrentGrid();
		rig.Render( 1 );
		const double cx = 0.5 * g.lx, cy = 0.5 * g.ly;
		//Outward flux of ball material through the circle, accumulated.
		const int bins = 72;
		std::vector< double > hist( bins, 0.0 );
		const int frames = 60;
		const Coils start = rig.plugin.CurrentCoils();
		auto coilAngle = [ & ]( const Coils& k ) { return std::atan2( k.y[ 0 ] - cy, k.x[ 0 ] - cx ); };
		double weighted = 0.0, weight = 0.0;
		for( int f = 0; f < frames; ++f )
		{
			rig.Render( 1 );
			const Snapshot s = Snapshot::Take( rig.plugin );
			double frameFlux = 0.0;
			for( int q = 0; q < 720; ++q )
			{
				const double a = 2.0 * kPi * ( q + 0.5 ) / 720.0;
				const double x = cx + rc * std::cos( a ), y = cy + rc * std::sin( a );
				const int i = static_cast< int >( x / s.dx ), j = static_cast< int >( y / s.dx );
				const double ur = s.U( i, j ) * std::cos( a ) + s.V( i, j ) * std::sin( a );
				const double carried = s.d[ s.At( i, j ) + 3 ];//rho chi
				if( ur > 0.0 )
				{
					hist[ static_cast< int >( a / ( 2.0 * kPi ) * bins ) % bins ] += carried * ur;
					frameFlux += carried * ur;
				}
			}
			//The coils' angle, weighted by what leaked while they stood there:
			//a smeared histogram's centre is where they were on average.
			weighted += frameFlux * AngleDiff( coilAngle( rig.plugin.CurrentCoils() ), coilAngle( start ) );
			weight += frameFlux;
		}
		const double turned = weight > 0.0 ? weighted / weight : 0.0;
		//Peaks: local maxima of the 3-bin smoothed histogram, largest N,
		//refined by a parabola through the three bins.
		std::vector< double > sm( bins );
		for( int b = 0; b < bins; ++b )
			sm[ b ] = hist[ ( b + bins - 1 ) % bins ] + hist[ b ] + hist[ ( b + 1 ) % bins ];
		std::vector< std::pair< double, double > > peaks;
		for( int b = 0; b < bins; ++b )
		{
			const double l = sm[ ( b + bins - 1 ) % bins ], m = sm[ b ], r = sm[ ( b + 1 ) % bins ];
			if( m > l && m >= r )
			{
				const double off = ( l - r ) / ( 2.0 * ( l - 2.0 * m + r ) );
				peaks.push_back( { m, 2.0 * kPi * ( b + 0.5 + off ) / bins } );
			}
		}
		std::sort( peaks.begin(), peaks.end(), []( auto& a, auto& b ) { return a.first > b.first; } );
		std::vector< double > predicted = CuspAngles( start, cx, cy, rc );
		for( double& a : predicted )
			a += turned;
		if( perturb.cuspAtCoils )
			for( double& a : predicted )
				a += kPi / c.poles;
		double worst = 0.0;
		bool enough  = static_cast< int >( peaks.size() ) >= c.poles && static_cast< int >( predicted.size() ) == c.poles;
		if( enough )
			for( double a : predicted )
			{
				double best = 1e9;
				for( int q = 0; q < c.poles; ++q )
					best = std::min( best, std::abs( AngleDiff( peaks[ q ].second, a ) ) );
				worst = std::max( worst, best );
			}
		//Half a bin, plus the angle two cells subtend at the circle.
		const double tol = kPi / bins + 2.0 * g.dx / rc;
		//And the N peaks must be leaks, not ripple: the weakest of them at
		//least half again the strongest bin halfway between two cusps. (An
		//eight-pole cusp's gaps are narrow and its contrast is 2.6.)
		double between = 0.0;
		for( double a : predicted )
		{
			const int b = static_cast< int >( std::fmod( a + kPi / c.poles + 2.0 * kPi, 2.0 * kPi ) / ( 2.0 * kPi ) * bins ) % bins;
			between     = std::max( between, sm[ b ] );
		}
		const double weakest = enough ? peaks[ c.poles - 1 ].first : 0.0;
		if( c.spin > 0.0 && enough )
		{
			//Following, not merely sitting: the coils as they STARTED must
			//miss the leaks by more than the bound.
			double stale = 0.0;
			for( double a : CuspAngles( start, cx, cy, rc ) )
			{
				double best = 1e9;
				for( int q = 0; q < c.poles; ++q )
					best = std::min( best, std::abs( AngleDiff( peaks[ q ].second, a ) ) );
				stale = std::max( stale, best );
			}
			Check( stale > kPi / bins + 2.0 * g.dx / rc,
			       fmt( "N = 6 turning: the unturned coils' cusps miss the leaks by %.2f deg (must exceed the bound)",
			            stale * 180.0 / kPi ) );
		}
		if( std::getenv( "CT_HIST" ) )
			for( int b = 0; b < bins; ++b )
				Say( "    %5.1f deg %10.4g\n", ( b + 0.5 ) * 360.0 / bins, sm[ b ] );
		Check( enough && worst < tol && weakest > 1.5 * between,
		       fmt( "N = %d%s: %zu cusp angles predicted, the %d strongest leaks within %.2f deg of them (bound %.2f); "
		            "weakest leak %.1fx the flow between cusps (bound 1.5)",
		            c.poles, c.spin > 0.0 ? fmt( ", Wall, coils turning %.1f rad/tau_A (turned %.1f deg, leak-weighted)",
		                                         c.spin, turned * 180.0 / kPi ).c_str() : "",
		            predicted.size(), c.poles, worst * 180.0 / kPi, tol * 180.0 / kPi, weakest / std::max( between, 1e-30 ) ) );
	}
	return g_failures;
}

//===========================================================================
// --frozen
//===========================================================================
int RunFrozen( const Perturb& perturb )
{
	Say( "\n=== frozen: Bz / rho is carried with the fluid (2.5-D, no in-plane field)\n" );
	Rig rig;
	rig.Init( 256, 256 );
	rig.Quiet();
	rig.Set( PT_DETAIL, DetailParam( 256 ) );
	rig.Set( PT_SPEED, 0.0f );
	rig.Set( PT_BOUNDARY, 1.0f );
	rig.Set( PT_POLES, PolesParam( 0 ) );
	rig.Set( PT_GUIDE_FIELD, GuideParam( 1.0 ) );
	rig.Set( PT_BALL_SIZE, 1.0f );//the drive's window follows the ball's size
	rig.Set( PT_DRIVE, 0.7f );
	rig.Set( PT_DRIVE_SCALE, 0.6f );
	rig.Render( 1 );
	const Grid& g = rig.plugin.CurrentGrid();
	StateBuilder sb( g.nx, g.ny, kGamma );
	auto pattern = [ & ]( double x, double y ) { return 1.0 + 0.5 * std::sin( 2.0 * kPi * x ) * std::sin( 2.0 * kPi * y ); };
	for( int j = 0; j < g.ny; ++j )
		for( int i = 0; i < g.nx; ++i )
		{
			const double x = ( i + 0.5 ) * g.dx, y = ( j + 0.5 ) * g.dx;
			const double rho = 1.0 + 0.3 * std::cos( 2.0 * kPi * x ) * std::cos( 4.0 * kPi * y );
			const double q   = pattern( x, y );
			const double bz  = rho * q;
			//Total pressure uniform, so nothing moves but what the drive moves.
			sb.Set( i, j, rho, 0, 0, 0, 3.0 - 0.5 * bz * bz, 0, 0, bz, q, 0, 0, 0 );
		}
	sb.Load( rig.plugin );
	rig.plugin.StepForTest( 1.5 );
	const Snapshot s = Snapshot::Take( rig.plugin );
	double sa = 0, sb2 = 0, saa = 0, sbb = 0, sab = 0, n = 0, moved = 0, spread = 0;
	for( int j = 0; j < s.ny; ++j )
		for( int i = 0; i < s.nx; ++i )
		{
			const double x = ( i + 0.5 ) * s.dx, y = ( j + 0.5 ) * s.dx;
			const double carried = s.Bz( i, j ) / s.Rho( i, j );
			const double tracer  = perturb.frozenStatic ? pattern( x, y ) : s.Tracer( i, j, 0 );
			sa += carried;
			sb2 += tracer;
			saa += carried * carried;
			sbb += tracer * tracer;
			sab += carried * tracer;
			n += 1.0;
			moved += std::pow( carried - pattern( x, y ), 2 );
			spread += std::pow( pattern( x, y ) - 1.0, 2 );
		}
	const double cov  = sab / n - ( sa / n ) * ( sb2 / n );
	const double corr = cov / std::sqrt( ( saa / n - sa * sa / n / n ) * ( sbb / n - sb2 * sb2 / n / n ) );
	const double displaced = std::sqrt( moved / spread );
	Check( displaced > 0.3, fmt( "the drive moved the pattern: Bz/rho differs from where it started by %.2f of its own "
	                             "spread (at least 0.3, or the check proves nothing)",
	                             displaced ) );
	//Two advections of the same mass flux -- Bz by the induction equation
	//through HLLD, the tracer upwinded -- that differ only in their
	//numerical diffusion at the pattern's scale (32 cells a wavelength).
	Check( corr > 0.99, fmt( "Bz/rho against the tracer carried with the mass%s: correlation %.5f (bound 0.99)",
	                         perturb.frozenStatic ? " (expected: the pattern where it STARTED)" : "", corr ) );
	return g_failures;
}

//===========================================================================
// --quench
//===========================================================================
namespace
{
/// Toro's exact Riemann solver for the Euler equations: the right-moving
/// shock's speed when a slab of (rhoL, pL) at rest meets (rhoR, pR) at rest.
double ExactShockSpeed( double rhoL, double pL, double rhoR, double pR, double gamma )
{
	const double cL = std::sqrt( gamma * pL / rhoL ), cR = std::sqrt( gamma * pR / rhoR );
	auto f = [ & ]( double p, double rho, double pk, double c, double& d ) {
		if( p > pk )
		{
			const double A = 2.0 / ( ( gamma + 1.0 ) * rho ), B = ( gamma - 1.0 ) / ( gamma + 1.0 ) * pk;
			const double q = std::sqrt( A / ( p + B ) );
			d              = q * ( 1.0 - ( p - pk ) / ( 2.0 * ( B + p ) ) );
			return ( p - pk ) * q;
		}
		const double r = p / pk;
		d = 1.0 / ( rho * c ) * std::pow( r, -( gamma + 1.0 ) / ( 2.0 * gamma ) );
		return 2.0 * c / ( gamma - 1.0 ) * ( std::pow( r, ( gamma - 1.0 ) / ( 2.0 * gamma ) ) - 1.0 );
	};
	double p = 0.5 * ( pL + pR );
	for( int it = 0; it < 200; ++it )
	{
		double dL, dR;
		const double F = f( p, rhoL, pL, cL, dL ) + f( p, rhoR, pR, cR, dR );
		const double next = std::max( p - F / ( dL + dR ), 1e-12 );
		if( std::abs( next - p ) < 1e-14 * p )
			break;
		p = next;
	}
	return cR * std::sqrt( ( gamma + 1.0 ) / ( 2.0 * gamma ) * p / pR + ( gamma - 1.0 ) / ( 2.0 * gamma ) );
}
} // namespace

int RunQuench( const Perturb& perturb )
{
	Say( "\n=== quench: the coils let go; the plasma free-expands\n" );
	//1. The expansion, in the plane: a slab of hot plasma with no field
	//left, into the tenuous gas around it. Its front can never outrun the
	//vacuum escape speed 2 c_s / ( gamma - 1 ); into a gas of density
	//rho_b it is a shock whose speed the exact Riemann solution gives, and
	//it approaches the escape speed as rho_b falls towards the floor.
	const double gamma = perturb.quenchGamma > 0.0 ? perturb.quenchGamma : kGamma;
	for( double rhoB : { 1e-2, 1e-3 } )
	{
		const double pB = 1e-5;
		TestModel model;
		model.boundary           = 3;
		model.backgroundDensity  = static_cast< float >( rhoB );
		model.backgroundPressure = static_cast< float >( pB );
		Rig rig;
		rig.plugin.SetModelForTest( model );
		rig.Init( 128, 128 );
		rig.Quiet();
		rig.Set( PT_DETAIL, DetailParam( 512 ) );
		rig.Set( PT_SPEED, 0.0f );
		rig.Set( PT_POLES, PolesParam( 0 ) );
		rig.Set( PT_GUIDE_FIELD, 0.0f );
		rig.Render( 1 );
		const Grid& g = rig.plugin.CurrentGrid();
		StateBuilder sb( g.nx, g.ny, kGamma );
		for( int j = 0; j < g.ny; ++j )
			for( int i = 0; i < g.nx; ++i )
			{
				const bool hot = ( i + 0.5 ) * g.dx < 0.25;
				sb.Set( i, j, hot ? 1.0 : rhoB, 0, 0, 0, hot ? 1.0 : pB, 0, 0, 0 );
			}
		sb.Load( rig.plugin );
		auto front = [ & ]() {
			const Snapshot s = Snapshot::Take( rig.plugin );
			const int j      = s.ny / 2;
			for( int i = s.nx - 1; i >= 0; --i )
				if( s.Rho( i, j ) > 2.5 * rhoB )
					return ( i + 0.5 ) * s.dx;
			return 0.0;
		};
		//From 0.08, once the shock has formed and its structure is steady, to
		//0.18.
		rig.plugin.StepForTest( 0.08 );
		const double x1 = front();
		rig.plugin.StepForTest( 0.1 );
		const double x2    = front();
		const double speed = ( x2 - x1 ) / 0.1;
		const double escape = 2.0 * std::sqrt( gamma * 1.0 / 1.0 ) / ( gamma - 1.0 );
		const double shock  = ExactShockSpeed( 1.0, 1.0, rhoB, pB, gamma );
		//The front is found to within a cell at each end: 2 dx over 0.1.
		const double tol = 2.0 * g.dx / 0.1;
		Check( speed <= escape + tol && std::abs( speed - shock ) < tol,
		       fmt( "into rho_b = %.0e: front at %.4f %s escape speed %.4f; exact Riemann shock %.4f (%.2f%% of escape; "
		            "measured error %.4f, bound 2 dx / dt = %.4f)",
		            rhoB, speed, speed <= escape + tol ? "<=" : "EXCEEDS", escape, shock, 100 * shock / escape,
		            std::abs( speed - shock ), tol ) );
	}

	//2. The quench itself: the default bottle, then Quench. The coils'
	//strength decays as exp( -t / kQuenchTime ) and the ball, let go, grows.
	{
		Rig rig;
		rig.Init( 320, 180 );
		rig.Set( PT_DRIVE, 0.0f );
		rig.Set( PT_CURVATURE, 0.0f );
		rig.Render( 120 );
		auto ballArea = [ & ]() {
			const Snapshot s = Snapshot::Take( rig.plugin );
			double a         = 0.0;
			for( int j = 0; j < s.ny; ++j )
				for( int i = 0; i < s.nx; ++i )
					if( s.Tracer( i, j, 3 ) > 0.5 )
						a += s.dx * s.dx;
			return a;
		};
		const double area0  = ballArea();
		const double guide0 = rig.plugin.CurrentCoils().guide;
		const double t0     = rig.plugin.SimTime();
		rig.Set( PT_QUENCH, 1.0f );
		rig.Render( 120 );
		const double t1       = rig.plugin.SimTime();
		const double expected = guide0 * std::exp( -( t1 - t0 ) / kQuenchTime );
		const double guide1   = rig.plugin.CurrentCoils().guide;
		const double area1    = ballArea();
		Check( std::abs( guide1 - expected ) < 1e-6 * guide0 + 1e-9,
		       fmt( "the guide field decays as exp( -t / %.2f ): %.6f after %.3f tau_A, expected %.6f", kQuenchTime,
		            guide1, t1 - t0, expected ) );
		Check( area1 > 1.5 * area0, fmt( "the ball let go: its area %.4f -> %.4f (x%.2f; at least x1.5)", area0, area1,
		                                 area1 / area0 ) );
	}
	return g_failures;
}

//===========================================================================
// --resist
//===========================================================================
int RunResist( const Perturb& perturb )
{
	Say( "\n=== resist: with eta > 0 the field diffuses into a static plasma on L^2 / eta\n" );
	const double eta = 1e-3, width = 0.08, bump = 0.05;
	Rig rig;
	rig.Init( 256, 256 );
	rig.Quiet();
	rig.Set( PT_DETAIL, DetailParam( 256 ) );
	rig.Set( PT_SPEED, 0.0f );
	//Wall. Open's margin is an absorbing layer that relaxes towards the
	//plugin's own ambient plasma (rho 0.1, p 0.005), and this plasma is
	//rho 100, p 1000: under Open the frame's edge becomes a Riemann problem
	//whose rarefaction reaches the centre long before the field has diffused
	//(the bump's peak read -0.996). The bump is 20 cells wide and 128 from
	//any wall, where Bz is the coils' own 1 to e^-(1.6/0.08)^2.
	rig.Set( PT_BOUNDARY, 1.0f );
	rig.Set( PT_POLES, PolesParam( 0 ) );
	rig.Set( PT_FIELD, FieldParam( 1.0 ) );
	rig.Set( PT_GUIDE_FIELD, GuideParam( 1.0 ) );
	rig.Set( PT_RESISTIVITY, ResistivityParam( eta ) );
	rig.Render( 1 );
	const Grid& g = rig.plugin.CurrentGrid();
	StateBuilder sb( g.nx, g.ny, kGamma );
	for( int j = 0; j < g.ny; ++j )
		for( int i = 0; i < g.nx; ++i )
		{
			const double x = ( i + 0.5 ) * g.dx - 0.5 * g.lx, y = ( j + 0.5 ) * g.dx - 0.5 * g.ly;
			const double bz = 1.0 + bump * std::exp( -( x * x + y * y ) / ( width * width ) );
			//Total pressure uniform: a static plasma. As B diffuses the
			//pressure hole that balanced it stays put, and the gas has to
			//squeeze by dp / ( gamma p ) to keep the total level -- which
			//compresses the frozen Bz by the same fraction. At p = 10 that
			//was 3% of the bump (the first version). At p = 1000 it is 1e-5;
			//rho = 100 keeps the sound speed, and the step, what it was.
			sb.Set( i, j, 100.0, 0, 0, 0, 1000.0 - 0.5 * bz * bz, 0, 0, bz );
		}
	sb.Load( rig.plugin );
	//The time for 4 eta t to equal w^2: the bump halves.
	const double t = width * width / ( 4.0 * eta );
	rig.plugin.StepForTest( t );
	const Snapshot s = Snapshot::Take( rig.plugin );
	const double peak     = s.Bz( g.nx / 2, g.ny / 2 ) - 1.0;//the cell nearest the centre
	const double off      = 0.5 * g.dx;
	const double w2       = width * width + 4.0 * perturb.resistFactor * eta * t;
	const double expected = bump * width * width / w2 * std::exp( -2.0 * off * off / w2 );
	//The five-point Laplacian's truncation on a Gaussian w/dx = 20 cells
	//wide is ( dx / w )^2 / 12 of it; the flows second order in the bump are
	//(bump)^2. 1% covers both with room.
	Check( std::abs( peak / expected - 1.0 ) < 0.01,
	       fmt( "after w^2 / ( 4 eta ) = %.2f tau_A: the bump's peak %.6f against the diffusion equation's %.6f%s (%.3f%%; bound 1%%)",
	            t, peak, expected, perturb.resistFactor != 1.0 ? " (with eta doubled)" : "",
	            100 * std::abs( peak / expected - 1.0 ) ) );
	return g_failures;
}

//===========================================================================
// --floors
//===========================================================================
int RunFloors( const Perturb& perturb )
{
	Say( "\n=== floors: how often they fire on the default look, and what happens without them\n" );
	{
		TestModel model;
		model.floors = !perturb.floorsOff;
		Rig rig;
		rig.plugin.SetModelForTest( model );
		rig.Init( 320, 180 );
		rig.Render( 600 );
		const Grid& g        = rig.plugin.CurrentGrid();
		const double cellSteps = static_cast< double >( g.nx ) * g.ny * rig.plugin.SubstepsTaken();
		const double fired   = rig.plugin.FloorHits() / cellSteps;
		const double entropy = rig.plugin.EntropyCells() / cellSteps;
		const bool finite    = Snapshot::Take( rig.plugin ).Finite();
		Check( finite && fired < 1e-4,
		       fmt( "600 frames of the default look: floors fired on %.2e of cell-steps (bound 1e-4); %.1f%% took the "
		            "dual-energy switch (the low-beta background)",
		            fired, 100 * entropy ) );
	}
	//Where the floors are needed: Einfeldt's 1-2-3 problem, two halves of a
	//plasma flying apart at Mach 30 each way. A gap opens once they separate
	//faster than 4 c / ( gamma - 1 ) = 6 c, and at Mach 30 it is a deep
	//vacuum; MUSCL-Hancock is not positivity-preserving there. (At Mach 3 and
	//5 the scheme's own diffusion kept the gap above 1e-3 and nothing was
	//needed: the negative control could not fail, so those were not tests.)
	{
		TestModel model;
		model.floors   = !perturb.floorsOff;
		model.boundary = 3;
		Rig rig;
		rig.plugin.SetModelForTest( model );
		rig.Init( 128, 128 );
		rig.Quiet();
		rig.Set( PT_DETAIL, DetailParam( 256 ) );
		rig.Set( PT_SPEED, 0.0f );
		rig.Set( PT_POLES, PolesParam( 0 ) );
		rig.Set( PT_GUIDE_FIELD, 0.0f );
		rig.Render( 1 );
		const Grid& g = rig.plugin.CurrentGrid();
		StateBuilder sb( g.nx, g.ny, kGamma );
		const double c = std::sqrt( kGamma * 0.4 );
		for( int j = 0; j < g.ny; ++j )
			for( int i = 0; i < g.nx; ++i )
				sb.Set( i, j, 1.0, ( i + 0.5 ) * g.dx < 0.5 ? -30.0 * c : 30.0 * c, 0, 0, 0.4, 0, 0, 0 );
		sb.Load( rig.plugin );
		const int steps  = rig.plugin.StepForTest( 0.1 );
		const Snapshot s = Snapshot::Take( rig.plugin );
		bool sane        = steps > 0 && s.Finite();
		double rmin = 1e30, pmin = 1e30;
		for( int j = 0; j < s.ny && sane; ++j )
			for( int i = 0; i < s.nx; ++i )
			{
				rmin = std::min( rmin, s.Rho( i, j ) );
				pmin = std::min( pmin, s.P( i, j ) );
			}
		sane = sane && rmin > 0.0 && pmin > -1e-6;
		Check( sane, fmt( "the 1-2-3 problem at Mach 30%s: %s (min rho %.2e, min p %.2e; floors fired %.0f times)",
		                  perturb.floorsOff ? ", FLOORS OFF" : "", sane ? "finite and positive" : "BLEW UP", rmin, pmin,
		                  rig.plugin.FloorHits() ) );
	}
	return g_failures;
}

//===========================================================================
// --mutation
//===========================================================================
int RunMutation( const Perturb& )
{
	//One character, in the shipped limiter: the monotonised-central slope's
	//centred difference, ( dl + dr ) / 2, becomes ( dl - dr ) / 2. Brio-Wu and
	//the Alfven wave must now fail, which proves the harness drives the shader
	//the plugin compiles and not a copy of it.
	const std::string find = "0.5 * abs( l + r )", replace = "0.5 * abs( l - r )";
	Say( "\n=== mutation: one character changed in the shipped shader -- \"%s\" -> \"%s\"\n", find.c_str(), replace.c_str() );
	int present = 0;
	for( const ProgramSource& source : AllSources() )
		if( source.fragment.find( find ) != std::string::npos )
			++present;
	Check( present > 0, fmt( "the text to mutate is in %d of the shipped programs", present ) );
	ContainmentPlugin::SetShaderMutationForTest( find, replace );
	const int before = g_failures;
	g_failures       = 0;
	RunBrioWu( Perturb() );
	const int briowu = g_failures;
	g_failures       = 0;
	RunAlfven( Perturb() );
	const int alfven = g_failures;
	g_failures       = before;
	ContainmentPlugin::SetShaderMutationForTest( "", "" );
	Check( briowu > 0 && alfven > 0, fmt( "against the mutated shader --briowu failed %d checks and --alfven %d (each must fail "
	                                      "at least one)", briowu, alfven ) );
	return g_failures;
}



//===========================================================================
// --open
//===========================================================================
namespace
{
struct LongRun
{
	double tau, floorFraction, rhoMin, bRatio, bRatio0, emissionPeak0, emissionPeak;
	bool finite;
};

LongRun RunLong( bool legacy, double tauTarget )
{
	TestModel model;
	model.legacyOpen = legacy;
	Rig rig;
	rig.plugin.SetModelForTest( model );
	rig.Init( 320, 180 );
	rig.Set( PT_BOUNDARY, static_cast< float >( Boundary::Open ) );
	rig.Set( PT_SPEED, ParamFromSpeed( 0.84f ) );//0.014 tau_A a frame
	LongRun out {};
	//|B| against the coils' own field over the same cells at the same moment.
	//The coils turn, and their field in the grid's corners -- nearest the
	//coil circle -- swings by 2x as a coil passes; a runaway is field the
	//coils are not making.
	auto measure = [ & ]( double& bratio, double& rmin, double& epeak ) {
		const Snapshot s   = Snapshot::Take( rig.plugin );
		const Coils& coils = rig.plugin.CurrentCoils();
		out.finite         = s.Finite();
		double bmax = 0.0, vmax = 0.0;
		rmin = 1e30;
		for( int j = 0; j < s.ny; ++j )
			for( int i = 0; i < s.nx; ++i )
			{
				bmax = std::max( bmax, std::sqrt( s.Bx( i, j ) * s.Bx( i, j ) + s.By( i, j ) * s.By( i, j )
				                                  + s.Bz( i, j ) * s.Bz( i, j ) ) );
				double bx, by, bz;
				VacuumField( coils, ( i + 0.5 ) * s.dx, ( j + 0.5 ) * s.dx, bx, by, bz );
				vmax = std::max( vmax, std::sqrt( bx * bx + by * by + bz * bz ) );
				rmin = std::min( rmin, s.Rho( i, j ) );
			}
		bratio         = bmax / vmax;
		const Floats e = ReadTexture( rig.plugin.EmissionTextureID(), rig.plugin.EmissionWidth(), rig.plugin.EmissionHeight() );
		epeak          = 0.0;
		for( size_t k = 3; k < e.size(); k += 4 )
			epeak = std::max( epeak, static_cast< double >( e[ k ] ) );
	};
	rig.Render( 60 );
	double r0 = 0.0;
	measure( out.bRatio0, r0, out.emissionPeak0 );
	out.bRatio = out.bRatio0;
	out.rhoMin = r0;
	while( rig.plugin.SimTime() < tauTarget && out.finite )
	{
		rig.Render( 100 );
		double b, r, e;
		measure( b, r, e );
		out.bRatio = std::max( out.bRatio, b );
		if( std::getenv( "CT_DBG" ) )
			Say( "    t %.2f  |B| max / the coils' %.3f  rho min %.4f  emission %.3f\n", rig.plugin.SimTime(), b, r, e );
		out.rhoMin       = std::min( out.rhoMin, r );
		out.emissionPeak = e;
	}
	const Grid& g     = rig.plugin.CurrentGrid();
	out.tau           = rig.plugin.SimTime();
	out.floorFraction = rig.plugin.FloorHits() / ( static_cast< double >( g.nx ) * g.ny * rig.plugin.SubstepsTaken() );
	return out;
}

/// What the boundary sends back, and whether it keeps it. A small blob of
/// compressed ambient plasma at the frame's centre, in a uniform axial field
/// (so the absorbing layer's reference IS the background), rings out as a
/// cylindrical fast wave. The same blob is run twice: in the frame as the
/// plugin bounds it, and in the same frame inside a plasma so big that
/// nothing from ITS edge can get back by the end. The runs are cell-for-cell
/// identical until the wave reaches the boundary, so their difference,
/// anywhere in the frame and at any time, is what the boundary sent back --
/// and nothing else: not the 2-D wake the wave leaves behind it, not the
/// scheme's dispersion.
///
/// (A plane pulse, the first version, is the wrong probe: it runs along the
/// top and bottom margins too, the layer damps it there as it should, and
/// the step that leaves diffracts into the frame. Nothing the ball sends out
/// runs along an edge from outside the frame.)
struct Echo
{
	double incident = 0.0;///< the wave's |dP| as it reaches the frame's edge
	double first    = 0.0;///< the largest echo in the frame once the wave has left it, over incident
	double left     = 0.0;///< what is still in the frame two crossings later, over incident
};

Echo Reflection( bool legacy )
{
	const double rho0 = kBackgroundDensity, p0 = kBackgroundPressure, b0 = 1.0, eps = 1e-3;
	const double cs2 = kGamma * p0 / rho0, cf = std::sqrt( ( kGamma * p0 + b0 * b0 ) / rho0 );
	const int cells = 256;
	const double dx = 1.0 / cells, sigma = 6.0 * dx;
	//The front clears the frame at its half-diagonal; an echo sent back from
	//anywhere on the edge has crossed the whole frame by tFirst; two more
	//crossings on, tLeft, an absorbing edge has had two more chances to take
	//what is left out.
	const double tClear = ( std::sqrt( 0.5 ) + 4.0 * sigma ) / cf;
	const double tFirst = tClear + std::sqrt( 2.0 ) / cf;
	const double tLeft  = tFirst + 2.0 * std::sqrt( 2.0 ) / cf;
	std::vector< double > times;
	for( int k = 1; k <= 4; ++k )
		times.push_back( 0.1 * k / cf );//the front near r = 0.4..0.5: the incident amplitude
	for( int k = 0; k <= 12; ++k )
		times.push_back( tClear + ( tFirst - tClear ) * k / 12.0 );
	times.push_back( tLeft );

	struct Run
	{
		std::vector< std::vector< double > > dP;//per time, frame cells
		int fx = 0, fy = 0;
	};
	auto run = [ & ]( const TestModel& model ) {
		Run out;
		Rig rig;
		rig.plugin.SetModelForTest( model );
		rig.Init( 256, 256 );
		rig.Quiet();
		rig.Set( PT_DETAIL, DetailParam( cells ) );
		rig.Set( PT_SPEED, 0.0f );
		rig.Set( PT_BOUNDARY, static_cast< float >( Boundary::Open ) );
		rig.Set( PT_POLES, PolesParam( 0 ) );
		rig.Set( PT_FIELD, FieldParam( 1.0 ) );
		rig.Set( PT_GUIDE_FIELD, GuideParam( 1.0 ) );
		rig.Render( 1 );
		const Grid& g = rig.plugin.CurrentGrid();
		out.fx = g.fx;
		out.fy = g.fy;
		const double cx = ( g.ox + 0.5 * g.fx ) * g.dx, cy = ( g.oy + 0.5 * g.fy ) * g.dx;
		StateBuilder sb( g.nx, g.ny, kGamma );
		for( int j = 0; j < g.ny; ++j )
			for( int i = 0; i < g.nx; ++i )
			{
				const double x = ( i + 0.5 ) * g.dx - cx, y = ( j + 0.5 ) * g.dx - cy;
				const double dr = eps * rho0 * std::exp( -0.5 * ( x * x + y * y ) / ( sigma * sigma ) );
				//An isentropic compression at rest, frozen flux: pure fast mode.
				sb.Set( i, j, rho0 + dr, 0, 0, 0, p0 + cs2 * dr, 0, 0, b0 * ( 1.0 + dr / rho0 ) );
			}
		sb.Load( rig.plugin );
		double now = 0.0;
		for( double t : times )
		{
			rig.plugin.StepForTest( t - now );
			now              = t;
			const Snapshot s = Snapshot::Take( rig.plugin );
			std::vector< double > frame;
			frame.reserve( static_cast< size_t >( g.fx ) * g.fy );
			for( int j = 0; j < g.fy; ++j )
				for( int i = 0; i < g.fx; ++i )
					frame.push_back( ( s.P( g.ox + i, g.oy + j ) - p0 )
					                 + 0.5 * ( s.Bz( g.ox + i, g.oy + j ) * s.Bz( g.ox + i, g.oy + j ) - b0 * b0 ) );
			out.dP.push_back( std::move( frame ) );
		}
		return out;
	};

	TestModel bounded;
	bounded.legacyOpen = legacy;
	TestModel unbounded;
	//A margin the wave cannot cross and come back through by tLeft, with no
	//layer in it: the frame inside an unbounded plasma.
	unbounded.marginFraction = static_cast< float >( 0.5 * ( cf * tLeft - 0.5 ) + 8.0 * dx );
	unbounded.spongeEFolds   = 0.0f;
	const Run a = run( bounded ), b = run( unbounded );
	Echo echo;
	if( a.fx != b.fx || a.fy != b.fy )
	{
		echo.first = echo.left = 1e30;
		return echo;
	}
	auto worstAt = [ & ]( size_t k ) {
		double worst = 0.0;
		for( size_t o = 0; o < a.dP[ k ].size(); ++o )
			worst = std::max( worst, std::abs( a.dP[ k ][ o ] - b.dP[ k ][ o ] ) );
		return worst;
	};
	//The incident amplitude: the unbounded run's largest |dP| within 0.1 of
	//the frame's inscribed circle, from the four snapshots with the front
	//between r = 0.4 and 0.5 (they are 0.1 apart in radius, so the front is
	//in that band in at least one).
	for( size_t k = 0; k < 4; ++k )
		for( int j = 0; j < a.fy; ++j )
			for( int i = 0; i < a.fx; ++i )
			{
				const double r = std::hypot( ( i + 0.5 - 0.5 * a.fx ) * dx, ( j + 0.5 - 0.5 * a.fy ) * dx );
				if( r > 0.4 && r <= 0.5 )
					echo.incident = std::max( echo.incident, std::abs( b.dP[ k ][ static_cast< size_t >( j ) * a.fx + i ] ) );
			}
	for( size_t k = 4; k + 1 < times.size(); ++k )
		echo.first = std::max( echo.first, worstAt( k ) );
	echo.left  = worstAt( times.size() - 1 );
	echo.first /= std::max( echo.incident, 1e-30 );
	echo.left /= std::max( echo.incident, 1e-30 );
	return echo;
}
} // namespace

int RunOpen( const Perturb& perturb )
{
	Say( "\n=== open: the frame is a window onto a bigger bottle -- waves leave, nothing drains, it lasts\n" );
	const bool legacy = perturb.legacyOpen;
	const Echo echo   = Reflection( legacy );
	Say( "  note  %sa cylindrical fast wave from the frame's centre: the first echo, once the wave has left the frame, "
	     "is %.2f%% of its amplitude at the edge (a low-frequency wave meeting the layer obliquely)\n",
	     legacy ? "OLD ZERO-GRADIENT GHOST, " : "", 100 * echo.first );
	Check( echo.left < echo.first * echo.first,
	       fmt( "%sthe echo leaves: two frame crossings later %.3f%% is still in the frame (bound: the first echo "
	            "squared, %.3f%%)",
	            legacy ? "OLD ZERO-GRADIENT GHOST, " : "", 100 * echo.left, 100 * echo.first * echo.first ) );
	const LongRun run = RunLong( legacy, 20.0 );
	//Sound: finite; floors as rare as `--floors` demands; the field no more
	//than twice the coils' own peak on the grid (a runaway goes past that in
	//a few tau_A); nowhere emptier than a tenth of the ambient density; and
	//the fuelled ball still there -- its emission at least a quarter of what
	//it was once it had settled.
	const bool ok = run.finite && run.tau >= 20.0 && run.floorFraction < 1e-4 && run.bRatio < 2.0
	                && run.rhoMin > 0.1 * kBackgroundDensity && run.emissionPeak > 0.25 * run.emissionPeak0;
	Check( ok, fmt( "%sthe default look for %.1f tau_A: finite %s, floors on %.1e of cell-steps (bound 1e-4), |B| at "
	                "most %.2f of the coils' own peak at the same moment (%.2f at the start; bound 2), rho min %.4f "
	                "(bound %.3f), the ball's emission %.3f of its settled %.3f (bound a quarter)",
	                legacy ? "OLD ZERO-GRADIENT GHOST, " : "", run.tau, run.finite ? "yes" : "NO", run.floorFraction,
	                run.bRatio, run.bRatio0, run.rhoMin, 0.1 * kBackgroundDensity, run.emissionPeak, run.emissionPeak0 ) );
	return g_failures;
}

//===========================================================================
// --presets
//===========================================================================
int RunPresets( const Perturb& )
{
	Say( "\n=== presets: each row sets exactly its values, Custom leaves the controls alone, row 1 is the defaults\n" );
	ContainmentPlugin plugin;
	const unsigned int targets[ presets::kParamCount ] = {
		PT_BALL_SIZE, PT_TEMPERATURE, PT_PROFILE, PT_FEED,    PT_CLIP_HEATS, PT_FUEL,        PT_FIELD,  PT_GUIDE_FIELD,
		PT_POLES,     PT_COIL_RADIUS, PT_COIL_SPIN, PT_CURVATURE, PT_QUENCH, PT_BOUNDARY,    PT_RESISTIVITY, PT_COOLING,
		PT_DRIVE,     PT_DRIVE_SCALE, PT_EXPOSURE, PT_TINT,    PT_RAMP,       PT_GLOW,        PT_FIELD_LINES, PT_LINE_COUNT,
	};
	const bool discrete[ presets::kParamCount ] = { false, false, true,  false, false, false, false, false,
		                                            true,  false, false, false, true,  true,  false, false,
		                                            false, false, false, false, true,  false, false, false };
	//The defaults, as the constructor set them.
	std::vector< float > defaults( PT_COUNT );
	for( unsigned int i = 0; i < PT_COUNT; ++i )
		defaults[ i ] = plugin.GetFloatParameter( i );
	int differ = 0;
	for( int c = 0; c < presets::kParamCount; ++c )
		if( defaults[ targets[ c ] ] != presets::kPresets[ 0 ].v[ c ] )
			++differ;
	Check( differ == 0, fmt( "row 1, \"%s\", is the constructor's defaults in all %d columns (%d differ)",
	                         presets::kPresets[ 0 ].name, presets::kParamCount, differ ) );
	//An operator's own values, all different from every row's.
	for( int c = 0; c < presets::kParamCount; ++c )
		plugin.SetFloatParameter( targets[ c ], discrete[ c ] ? 0.0f : 0.123f );
	plugin.SetFloatParameter( PT_BALL_X, 0.321f );
	for( int r = 0; r <= presets::kCount; ++r )
	{
		plugin.SetFloatParameter( PT_PRESET, static_cast< float >( r ) );
		int wrong = 0, fractions = 0;
		for( int c = 0; c < presets::kParamCount; ++c )
		{
			const float want = r == 0 ? ( discrete[ c ] ? 0.0f : 0.123f ) : presets::kPresets[ r - 1 ].v[ c ];
			if( plugin.EffectiveForTest( targets[ c ] ) != want )
				++wrong;
			if( r > 0 && discrete[ c ] && want != std::round( want ) )
				++fractions;
		}
		//What a row does not own stays the operator's.
		if( plugin.EffectiveForTest( PT_BALL_X ) != 0.321f )
			++wrong;
		Check( wrong == 0 && fractions == 0,
		       fmt( "%-20s every column as the row says, Ball X still the operator's (%d wrong, %d fractions in discrete "
		            "columns)",
		            r == 0 ? "Custom:" : ( std::string( presets::kPresets[ r - 1 ].name ) + ":" ).c_str(), wrong, fractions ) );
	}
	//And the host sees the rows by name.
	Check( std::string( plugin.GetParamName( PT_PRESET ) ) == "Preset",
	       fmt( "the dropdown is called \"%s\" and has %d rows plus Custom", plugin.GetParamName( PT_PRESET ), presets::kCount ) );
	return g_failures;
}

//===========================================================================
// The checks that need no GL context (--offline). A GitHub macOS runner cannot
// make an accelerated one, so these are what CI can actually run.
//===========================================================================

// --reference: the Brio-Wu reference the GPU is judged against, on its own.
int RunReference( const Perturb& perturb )
{
	const double gamma = 2.0, bn = 0.75, tEnd = 0.1;
	const double gammaRun = perturb.referenceGamma > 0.0 ? perturb.referenceGamma : gamma;
	Say( "\n=== reference: the double-precision Brio-Wu solver's fast heads against the closed-form fast speeds%s\n",
	     perturb.referenceGamma > 0.0 ? " (THE REFERENCE RUN WITH THE WRONG GAMMA)" : "" );
	BrioWuReference ref;
	ref.n     = 8192;
	ref.gamma = gammaRun;
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
	const Features R   = Locate( rx, rr, rp, rb );
	const double headL = 0.5 - tEnd * FastSpeed( 1.0, 1.0, bn, 1.0, 0.0, gamma );
	const double headR = 0.5 + tEnd * FastSpeed( 0.125, 0.1, bn, -1.0, 0.0, gamma );
	//The same bound as --briowu's: 4 cells for where a head is found at 0.1%
	//of the jump, plus 0.002. A rarefaction's head is a kink, where minmod
	//drops to first order and smears it diffusively over the run: a width of
	//sqrt( c dx t / 2 ) = 0.0033 at 8192 cells (c = 1.8, t = 0.1). The kink's
	//0.1% point moves by less than that width.
	Check( std::abs( R.fastLeft - headL ) < 4.0 / ref.n + 0.002 && std::abs( R.fastRight - headR ) < 4.0 / ref.n + 0.002,
	       fmt( "the reference's fast heads sit where the fast speeds put them (%.4f vs %.4f, %.4f vs %.4f; bound %.4f)",
	            R.fastLeft, headL, R.fastRight, headR, 4.0 / ref.n + 0.002 ) );
	return g_failures;
}

// --vacuum: the coils' closed-form field, in double.
int RunVacuum( const Perturb& perturb )
{
	Say( "\n=== vacuum: the coils' field is a vacuum field (no divergence, no current), a multipole, normalised%s\n",
	     perturb.vacuumSameSign ? " (EVERY COIL CARRYING THE SAME CURRENT)" : "" );
	const double cx = 0.9, cy = 0.5, radius = 1.2, field = 0.7, h = 1e-3;
	for( int poles : { 4, 6, 8, 12 } )
	{
		Coils coils = MakeCoils( poles, radius, 0.3, field, 0.0, 1.0, cx, cy );
		if( perturb.vacuumSameSign )
			for( int k = 0; k < coils.count; ++k )
				coils.current[ k ] = std::abs( coils.current[ k ] );
		auto b = [ & ]( double x, double y, double& bx, double& by ) {
			double bz;
			VacuumField( coils, x, y, bx, by, bz );
		};
		//div B and the in-plane curl (J_z), by central differences at points
		//inside the frame's reach, against what the differencing itself must
		//leave. Its truncation is h^2 / 6 times a third derivative somewhere
		//within h, and a line current's 1/d field has third derivatives up to
		//6 |I| / d^4: so at most |I| h^2 / ( d - h )^4 per derivative, two
		//derivatives, summed over the coils. Its rounding: each field value is a sum of N terms of up to
		//|I| / d, good to N u of that, divided by 2h, twice.
		Pcg random( 7 + poles );
		double worst = 0.0;
		const double u = std::numeric_limits< double >::epsilon();
		for( int n = 0; n < 400; ++n )
		{
			const double a = 2.0 * kPi * random.Uniform(), r = 0.95 * random.Uniform();
			const double x = cx + r * std::cos( a ), y = cy + r * std::sin( a );
			double bxp, byp, bxm, bym, bxu, byu, bxd, byd;
			b( x + h, y, bxp, byp );
			b( x - h, y, bxm, bym );
			b( x, y + h, bxu, byu );
			b( x, y - h, bxd, byd );
			double truncation = 0.0, scale = 0.0;
			for( int k = 0; k < coils.count; ++k )
			{
				const double d = std::hypot( x - coils.x[ k ], y - coils.y[ k ] );
				const double near = d - h;
				truncation += 2.0 * std::abs( coils.current[ k ] ) * h * h / ( near * near * near * near );
				scale += std::abs( coils.current[ k ] ) / d;
			}
			const double allowed = truncation + 2.0 * coils.count * u * scale / h;
			const double div     = ( bxp - bxm + byu - byd ) / ( 2.0 * h );
			const double cur     = ( byp - bym - bxu + bxd ) / ( 2.0 * h );
			worst = std::max( worst, std::max( std::abs( div ), std::abs( cur ) ) / allowed );
		}
		const double bound = 1.0;
		//Near the centre a 2m-pole field grows as r^( m - 1 ), m = N / 2. The
		//next term of the coils' expansion is ( r / R )^N of the first.
		const double gap = 0.3 + kPi / poles;
		double b1x, b1y, b2x, b2y;
		b( cx + 0.05 * std::cos( gap ), cy + 0.05 * std::sin( gap ), b1x, b1y );
		b( cx + 0.1 * std::cos( gap ), cy + 0.1 * std::sin( gap ), b2x, b2y );
		const double slope  = std::log( std::hypot( b2x, b2y ) / std::hypot( b1x, b1y ) ) / std::log( 2.0 );
		const double expect = poles / 2 - 1;
		const double sBound = 4.0 * poles * std::pow( 0.1 / radius, poles ) / std::log( 2.0 ) + 1e-9;
		//Normalised: |B| at radius 0.5 in the first gap is Field.
		double bnx, bny;
		b( cx + 0.5 * std::cos( gap ), cy + 0.5 * std::sin( gap ), bnx, bny );
		const double norm = std::hypot( bnx, bny ) / field - 1.0;
		const std::vector< double > cusps = CuspAngles( coils, cx, cy, 0.5 );
		Check( worst < bound && std::abs( slope - expect ) < sBound && std::abs( norm ) < 1e-12
		           && static_cast< int >( cusps.size() ) == poles,
		       fmt( "N = %d: div B and J_z at most %.3f of the difference's truncation (bound 1); |B| grows as r^%.6f "
		            "(r^%.0f expected, bound %.1e); |B| at 0.5 in the gap is Field to %.1e (bound 1e-12); %zu cusps",
		            poles, worst, slope, expect, sBound, std::abs( norm ), cusps.size() ) );
	}
	return g_failures;
}

// --names: what an operator and `--set` see.
int RunNames( const Perturb& perturb )
{
	Say( "\n=== names: every parameter has a unique, human name, and the display name fits the host's field\n" );
	ContainmentPlugin plugin;
	std::vector< std::string > names;
	int empty = 0, duplicates = 0;
	for( unsigned int i = 0; i < PT_COUNT; ++i )
	{
		const char* name = plugin.GetParamName( i );
		std::string n    = name ? name : "";
		if( perturb.namesDuplicate && i == PT_FEED )
			n = plugin.GetParamName( PT_FUEL );
		if( n.empty() )
			++empty;
		if( std::find( names.begin(), names.end(), n ) != names.end() )
			++duplicates;
		names.push_back( n );
	}
	Check( empty == 0 && duplicates == 0,
	       fmt( "%u parameters: %d unnamed, %d duplicate names%s", static_cast< unsigned int >( PT_COUNT ), empty,
	            duplicates, perturb.namesDuplicate ? " (Feed RENAMED to Fuel's name)" : "" ) );
	//The FFGL name field is char[ 16 ] and not null-terminated: a longer name
	//is truncated by the host without a word. The name itself lives in
	//PluginEntry.cpp, which only the bundle links; oxbow reads it back there.
	const std::string display = kDisplayName;
	Check( display.size() <= 16 && display.rfind( "SW ", 0 ) == 0,
	       fmt( "the display name \"%s\" is %zu bytes (bound 16) and carries the fleet's \"SW \" prefix", display.c_str(),
	            display.size() ) );
	return g_failures;
}

//===========================================================================
// --negative
//===========================================================================
int RunNegative( bool offline )
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
	add( "balance", RunBalance, "expect Bz inside = sqrt( B0^2 - p_in ), without the 2", []( Perturb& p ) { p.balanceNoTwo = true; } );
	add( "balance", RunBalance, "run with resistivity 1e-3: the field diffuses into the ball", []( Perturb& p ) { p.balanceEta = 1e-3; } );
	add( "rt", RunRT, "the wrong sign on Curvature", []( Perturb& p ) { p.curvatureSign = -1.0; } );
	add( "cusp", RunCusp, "expect the leaks towards the coils", []( Perturb& p ) { p.cuspAtCoils = true; } );
	add( "frozen", RunFrozen, "compare with where the pattern started", []( Perturb& p ) { p.frozenStatic = true; } );
	add( "quench", RunQuench, "predict with gamma 7/5", []( Perturb& p ) { p.quenchGamma = 1.4; } );
	add( "resist", RunResist, "expect eta twice as big", []( Perturb& p ) { p.resistFactor = 2.0; } );
	add( "floors", RunFloors, "run with the floors off", []( Perturb& p ) { p.floorsOff = true; } );
	add( "reference", RunReference, "run the reference with gamma 5/3 against gamma 2's fast speeds",
	     []( Perturb& p ) { p.referenceGamma = 5.0 / 3.0; } );
	add( "vacuum", RunVacuum, "every coil carrying the same current", []( Perturb& p ) { p.vacuumSameSign = true; } );
	add( "names", RunNames, "Feed given Fuel's name", []( Perturb& p ) { p.namesDuplicate = true; } );
	add( "open", RunOpen, "Open as it was: the zero-gradient ghost, no margin, no absorbing layer",
	     []( Perturb& p ) { p.legacyOpen = true; } );
	//The ones --offline may run: they make no GL context.
	if( offline )
		cases.erase( std::remove_if( cases.begin(), cases.end(),
		                             []( const Case& c ) {
			                             const std::string n = c.name;
			                             return n != "reference" && n != "vacuum" && n != "names";
		                             } ),
		             cases.end() );
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
