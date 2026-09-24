#pragma once

/**
    cttest's shared plumbing: reporting, the GL context, the rig that drives the
    REAL plugin class through ProcessOpenGL on a synthetic 60 fps clock, a
    snapshot of the plasma's state in double, and the PNG writer.
*/

#include "Containment.h"
#include "Controls.h"
#include "Physics.h"

#include <OpenGL/OpenGL.h>
#include <OpenGL/gl3.h>

#include <cmath>
#include <string>
#include <vector>

namespace cttest
{
using namespace containment;

constexpr double kPi = 3.14159265358979323846;
using Floats = std::vector< float >;

//---------------------------------------------------------------------------
// Reporting. Every check prints the numbers it compared.
//---------------------------------------------------------------------------
extern int g_failures;

/// `--size WxH` given with a check: every rig the check makes renders at this
/// raster instead of its own (0 = each check's own). The grid is Detail cells
/// on its short side whatever the raster; the raster sets the grid's aspect
/// and everything drawn. verify.sh runs every check twice: at each check's
/// own raster and at 320x180, CI's.
extern int g_rasterW, g_rasterH;
void ChooseRaster( int& width, int& height );
std::string fmt( const char* format, ... );
void Check( bool condition, const std::string& message );
void Say( const char* format, ... );
int Verdict();

/**
    Deliberate errors, for --negative. Each field, when set, makes one check
    either score the plugin against a model that is wrong by an amount the
    check should see, or run the plugin as a wrong model. A default one is the
    truth.
*/
struct Perturb
{
	bool hll              = false;///< run the plugin with HLL where HLLD is claimed
	bool glmOff           = false;///< run the plugin with no divergence cleaning
	bool floorsOff        = false;///< run the plugin with no floors
	double curvatureSign  = 1.0;  ///< --rt expects g with this sign
	bool alfvenOverRho    = false;///< --alfven expects B / rho, not B / sqrt( rho )
	bool balanceNoTwo     = false;///< --balance expects Bz_in = sqrt( B0^2 - p_in )
	double balanceEta     = 0.0;  ///< --balance runs with this resistivity (0: ideal, as claimed)
	double referenceGamma = 0.0;  ///< --reference runs its solver with this gamma (0: the true 2)
	bool vacuumSameSign   = false;///< --vacuum gives every coil the same current
	bool namesDuplicate   = false;///< --names sees Feed renamed to Fuel's name
	bool cuspAtCoils      = false;///< --cusp expects the leaks towards the coils
	bool conserveCooling  = false;///< --conserve runs with a real energy sink on
	bool frozenStatic     = false;///< --frozen compares with where the pattern STARTED
	double quenchGamma    = 0.0;  ///< --quench predicts with this gamma (0: the true one)
	double resistFactor   = 1.0;  ///< --resist expects this times eta
	bool stillUlp         = false;///< --still expects one ulp more in one pixel
	bool additiveGlow     = false;///< run the plugin's glow as an additive bloom
	bool legacyOpen       = false;///< run Open as it was: zero-gradient ghost, no absorbing layer
};

//---------------------------------------------------------------------------
// GL.
//---------------------------------------------------------------------------
CGLContextObj CreateContext();
GLuint MakeTexture( int width, int height, const float* pixels );
Floats ReadTexture( GLuint texture, int width, int height, GLenum format = GL_RGBA, int channels = 4 );

/// The card: a picture with colour and detail at every scale, so the clip has
/// something to be made of. Rows v = 0 first, the way GL stores a texture.
Floats BuildCard( int width, int height );

bool WritePng( const std::string& path, int width, int height, const Floats& rgba );

//---------------------------------------------------------------------------
// The audio the harness feeds (rosette's): a beat every half second.
//---------------------------------------------------------------------------
enum class AudioFeed
{
	Silence,
	Pulses
};

//---------------------------------------------------------------------------
// The rig: the real plugin, a float input texture and a float output
// framebuffer, on a synthetic 60 fps clock.
//---------------------------------------------------------------------------
struct Rig
{
	ContainmentPlugin plugin;
	int width = 0, height = 0;
	GLuint sourceTexture = 0, outputTexture = 0, outputFBO = 0;
	int frame            = 0;
	double fps           = 60.0;
	AudioFeed feed       = AudioFeed::Silence;

	ProcessOpenGLStruct process    = {};
	FFGLTextureStruct inputStruct  = {};
	FFGLTextureStruct* inputs[ 1 ] = { nullptr };

	~Rig();
	bool Init( int w, int h, const Floats* picture = nullptr );
	void Upload( const Floats& picture );
	void Set( unsigned int id, float value );
	void Press( unsigned int id );
	bool Render( int frames = 1 );
	Floats Output() const;

	/// A physics rig: the plasma alone. No feed, no drive, no spin, no
	/// curvature, no audio: the only thing happening is what the check sets.
	void Quiet();
};

/// The inverse mappings, so a check can ask for physical values by name.
float ParamGeometric( double value, double low, double high );
float DetailParam( int cells );
float PolesParam( int poles );
float CurvatureParam( double g );
float FieldParam( double field );
float GuideParam( double multiple );
float TemperatureParam( double beta );
float ResistivityParam( double eta );
float BallSizeParam( double radius );

//---------------------------------------------------------------------------
// The state, read back, in double.
//---------------------------------------------------------------------------
struct Snapshot
{
	int nx = 0, ny = 0;
	double dx = 0.0;
	double gamma = 5.0 / 3.0;
	Floats a, b, c, d;

	static Snapshot Take( const ContainmentPlugin& plugin, double gamma = 5.0 / 3.0 );

	size_t At( int i, int j ) const
	{
		return ( static_cast< size_t >( j ) * nx + i ) * 4;
	}
	double Rho( int i, int j ) const
	{
		return a[ At( i, j ) ];
	}
	double U( int i, int j ) const
	{
		return a[ At( i, j ) + 1 ] / a[ At( i, j ) ];
	}
	double V( int i, int j ) const
	{
		return a[ At( i, j ) + 2 ] / a[ At( i, j ) ];
	}
	double W( int i, int j ) const
	{
		return a[ At( i, j ) + 3 ] / a[ At( i, j ) ];
	}
	double Energy( int i, int j ) const
	{
		return b[ At( i, j ) ];
	}
	double Bx( int i, int j ) const
	{
		return b[ At( i, j ) + 1 ];
	}
	double By( int i, int j ) const
	{
		return b[ At( i, j ) + 2 ];
	}
	double Bz( int i, int j ) const
	{
		return b[ At( i, j ) + 3 ];
	}
	double Psi( int i, int j ) const
	{
		return c[ At( i, j ) ];
	}
	double P( int i, int j ) const;
	/// The tracer channel k (0..2) and the ball marker (3), per unit mass.
	double Tracer( int i, int j, int k ) const
	{
		return d[ At( i, j ) + k ] / a[ At( i, j ) ];
	}
	bool Finite() const;
	double Sum( int texture, int channel ) const;
};

/// Build a state from primitive values, one cell at a time.
struct StateBuilder
{
	int nx, ny;
	double gamma;
	Floats a, b, c, d;
	StateBuilder( int nx, int ny, double gamma );
	void Set( int i, int j, double rho, double u, double v, double w, double p, double bx, double by, double bz,
	          double t0 = 0.0, double t1 = 0.0, double t2 = 0.0, double t3 = 0.0 );
	void Load( ContainmentPlugin& plugin ) const
	{
		plugin.LoadStateForTest( a, b, c, d );
	}
};

} // namespace cttest
