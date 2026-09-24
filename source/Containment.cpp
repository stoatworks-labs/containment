#include "Containment.h"

#include "Diag.h"
#include "GLState.h"
#include "Presets.h"

//FFGLSDK.h includes every other scoped binding and omits this one (SDK
//b1afaf9), so it has to be asked for by name.
#include <ffglex/FFGLScopedFBOBinding.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>

using namespace ffglex;

namespace containment
{
namespace
{
constexpr double kPi = 3.14159265358979323846;

std::string glStringOrUnknown( GLenum name )
{
	const GLubyte* value = glGetString( name );
	return value ? reinterpret_cast< const char* >( value ) : "unknown";
}

const char* const kProfileNames[]  = { "Gaussian", "Top Hat" };
const char* const kPoleNames[]     = { "0", "4", "6", "8", "12" };
const char* const kBoundaryNames[] = { "Open", "Wall" };
const char* const kDetailNames[]   = { "128", "256", "512", "1024" };
const char* const kRampNames[]     = { "Hot", "Aurora" };
/// Texels per side one reduction pass folds into one. Small on purpose: a
/// fragment that loops over hundreds of texels runs them one after another,
/// and at the default grid the reduction is latency, not bandwidth. (16x16
/// blocks and a clock that finished the reduction in one fragment made a
/// substep 2.4x SLOWER than a chain of 4x4 passes.)
constexpr int kReduceBlock = 4;

const char* const kViewNames[]     = { "Picture", "Density", "Pressure", "|B|", "Beta", "Speed", "Div B" };

constexpr int kClockVotes = 4;

/// Seconds of host time one frame may advance by. The host's clock jumps
/// when a composition is scrubbed or the machine sleeps; a plasma that runs
/// a minute in one frame has simply lost everything it was doing.
constexpr double kMaxFrameDelta = 0.1;

std::string g_mutationFind;
std::string g_mutationReplace;

double wallSeconds()
{
	using namespace std::chrono;
	static const steady_clock::time_point start = steady_clock::now();
	return duration_cast< duration< double > >( steady_clock::now() - start ).count();
}

int optionIndex( float value, int count )
{
	return std::clamp( static_cast< int >( std::lround( value ) ), 0, count - 1 );
}

GLint location( GLuint program, const char* name )
{
	return glGetUniformLocation( program, name );
}
void uniform1i( GLuint p, const char* n, int v )
{
	glUniform1i( location( p, n ), v );
}
void uniform1f( GLuint p, const char* n, float v )
{
	glUniform1f( location( p, n ), v );
}
void uniform2f( GLuint p, const char* n, float a, float b )
{
	glUniform2f( location( p, n ), a, b );
}
void uniform3f( GLuint p, const char* n, float a, float b, float c )
{
	glUniform3f( location( p, n ), a, b, c );
}
void uniform2i( GLuint p, const char* n, int a, int b )
{
	glUniform2i( location( p, n ), a, b );
}

/// Bind textures to units 0..n-1, and unbind them all on the way out,
/// leaving unit 0 active. By hand rather than with the SDK's scoped bindings:
/// those clear to 0 on whichever unit is active when they exit, and with more
/// than two interleaved pairs that leaks a texture into the host's context
/// (millpond's trap). This plugin binds up to nine.
class BoundTextures
{
public:
	BoundTextures( std::initializer_list< GLuint > ids ) : count( static_cast< int >( ids.size() ) )
	{
		int unit = 0;
		for( GLuint id : ids )
		{
			glActiveTexture( static_cast< GLenum >( GL_TEXTURE0 + unit++ ) );
			glBindTexture( GL_TEXTURE_2D, id );
		}
		glActiveTexture( GL_TEXTURE0 );
	}
	~BoundTextures()
	{
		for( int unit = count - 1; unit >= 0; --unit )
		{
			glActiveTexture( static_cast< GLenum >( GL_TEXTURE0 + unit ) );
			glBindTexture( GL_TEXTURE_2D, 0 );
		}
		glActiveTexture( GL_TEXTURE0 );
	}
	BoundTextures( const BoundTextures& ) = delete;
	BoundTextures& operator=( const BoundTextures& ) = delete;

private:
	int count;
};

/// Samplers named in order onto units 0..n-1.
void samplers( GLuint program, std::initializer_list< const char* > names )
{
	int unit = 0;
	for( const char* name : names )
		uniform1i( program, name, unit++ );
}

/// A 1-D Gaussian of this sigma, in cells: tap count and the weights for
/// |offset| = 0..taps, normalised over the full, symmetric kernel.
int gaussianWeights( double sigma, std::vector< float >& weights )
{
	const int taps = std::clamp( static_cast< int >( std::ceil( 3.0 * sigma ) ), 1, kMaxBlurTaps );
	weights.assign( static_cast< size_t >( kMaxBlurTaps + 1 ), 0.0f );
	double sum = 0.0;
	std::vector< double > w( static_cast< size_t >( taps + 1 ) );
	for( int k = 0; k <= taps; ++k )
	{
		w[ k ] = std::exp( -0.5 * k * k / std::max( sigma * sigma, 1e-12 ) );
		sum += k == 0 ? w[ k ] : 2.0 * w[ k ];
	}
	for( int k = 0; k <= taps; ++k )
		weights[ k ] = static_cast< float >( w[ k ] / sum );
	return taps;
}
} // namespace

/// Which ParamId each preset column drives.
constexpr ParamId kPresetTarget[ presets::kParamCount ] = {
	PT_BALL_SIZE, PT_TEMPERATURE, PT_PROFILE, PT_FEED,    PT_CLIP_HEATS, PT_FUEL,        PT_FIELD,  PT_GUIDE_FIELD,
	PT_POLES,     PT_COIL_RADIUS, PT_COIL_SPIN, PT_CURVATURE, PT_QUENCH, PT_BOUNDARY,    PT_RESISTIVITY, PT_COOLING,
	PT_DRIVE,     PT_DRIVE_SCALE, PT_EXPOSURE, PT_TINT,    PT_RAMP,       PT_GLOW,        PT_FIELD_LINES, PT_LINE_COUNT,
};

float ContainmentPlugin::P( unsigned int index ) const
{
	const int preset = std::clamp( static_cast< int >( std::lround( params[ PT_PRESET ] ) ), 0, presets::kCount );
	if( preset > 0 && index != PT_PRESET )
		for( int c = 0; c < presets::kParamCount; ++c )
			if( kPresetTarget[ c ] == index )
				return presets::kPresets[ preset - 1 ].v[ c ];
	return index < PT_COUNT ? params[ index ] : 0.0f;
}

static_assert( PT_COUNT - PT_ABOUT_TEXT == stoatworks::about::kParamCount,
               "the About run no longer matches StoatworksAbout.h -- add or remove a PT_ABOUT_BUTTON_n to match" );

//---------------------------------------------------------------------------
void ContainmentPlugin::SetShaderMutationForTest( const std::string& find, const std::string& replace )
{
	g_mutationFind    = find;
	g_mutationReplace = replace;
}

//---------------------------------------------------------------------------
ContainmentPlugin::ContainmentPlugin()
{
	SetMinInputs( 1 );
	SetMaxInputs( 1 );
	SetTimeSupported( true );

	//-------------------------------------------------------------------
	// Defaults: a ball held mostly by the guide field and shaped by a weak
	// six-pole cusp, in a conducting vessel (Wall: the field the ball
	// compresses stays in the frame rather than leaking out of the nearest
	// edge), the coils turning slowly, the curvature and the stirring set
	// just past where the edge stays smooth -- so dropping the effect on a
	// layer shows a live, writhing plasma straight away. The clip keeps its
	// own colours (Temperature Tint 0); docs/green-orb.preset is the green look.
	//-------------------------------------------------------------------
	params[ PT_PRESET ]      = 0.0f;//Custom: the controls are the truth
	params[ PT_IGNITE ]      = 0.0f;
	params[ PT_BALL_SIZE ]   = 0.389f;//0.18 frame heights
	params[ PT_BALL_X ]      = 0.5f;
	params[ PT_BALL_Y ]      = 0.5f;
	params[ PT_TEMPERATURE ] = 0.315f;//beta0 = 0.7
	params[ PT_PROFILE ]     = static_cast< float >( Profile::Gaussian );
	params[ PT_PELLET ]      = 0.0f;
	params[ PT_FEED ]        = 0.1f;
	params[ PT_CLIP_HEATS ]  = 0.0f;
	params[ PT_FUEL ]        = 0.15f; //0.3 per tau_A: an open bottle that does not run down

	params[ PT_FIELD ]       = 0.25f; //0.5 B_ref at the rim
	params[ PT_GUIDE_FIELD ] = 0.7f;  //Bz = 1.4 Field: the guide field holds the ball, the cusp shapes it
	params[ PT_POLES ]       = 2.0f;  //six
	params[ PT_COIL_RADIUS ] = 0.105f;//1.3 frame heights
	params[ PT_COIL_SPIN ]   = 0.55f; //0.1 rad per tau_A
	params[ PT_CURVATURE ]   = 0.316f;//g_eff 0.3
	params[ PT_QUENCH ]      = 0.0f;
	params[ PT_BOUNDARY ]    = static_cast< float >( Boundary::Open );

	params[ PT_SPEED ]       = ParamFromSpeed( 0.3f );
	params[ PT_RESISTIVITY ] = 0.0f;
	params[ PT_COOLING ]     = 0.0f;
	params[ PT_DETAIL ]      = 1.0f;  //256 cells on the short side
	params[ PT_DRIVE ]       = 0.274f;//0.3 rms
	params[ PT_DRIVE_SCALE ] = 0.442f;//0.15 frame heights

	params[ PT_AUDIO_HEAT ]    = 0.0f;
	params[ PT_AUDIO_PELLETS ] = 0.0f;

	params[ PT_EXPOSURE ]    = 0.4167f;//+1 EV
	params[ PT_TINT ]        = 0.0f;
	params[ PT_RAMP ]        = static_cast< float >( Ramp::Hot );
	params[ PT_GLOW ]        = 0.333f; //0.3 of the light into the glare
	params[ PT_FIELD_LINES ] = 0.0f;
	params[ PT_LINE_COUNT ]  = 0.273f; //16
	params[ PT_VIEW ]        = static_cast< float >( View::Picture );
	params[ PT_MIX ]         = 1.0f;

	//-------------------------------------------------------------------
	// Declaration. Every numeric parameter is 0..1; Controls.cpp maps them.
	//-------------------------------------------------------------------
	{
		SetOptionParamInfo( PT_PRESET, "Preset", presets::kCount + 1, 0.0f );
		SetParamElementInfo( PT_PRESET, 0, "Custom", 0.0f );
		for( int i = 0; i < presets::kCount; ++i )
			SetParamElementInfo( PT_PRESET, static_cast< unsigned int >( i + 1 ), presets::kPresets[ i ].name,
			                     static_cast< float >( i + 1 ) );
	}
	SetParamInfo( PT_IGNITE, "Ignite", FF_TYPE_EVENT, false );
	SetParamInfof( PT_BALL_SIZE, "Ball Size", FF_TYPE_STANDARD );
	SetParamInfof( PT_BALL_X, "Ball X", FF_TYPE_XPOS );
	SetParamInfof( PT_BALL_Y, "Ball Y", FF_TYPE_YPOS );
	SetParamInfof( PT_TEMPERATURE, "Temperature", FF_TYPE_STANDARD );
	SetOptionParamInfo( PT_PROFILE, "Profile", static_cast< int >( Profile::Count ), params[ PT_PROFILE ] );
	for( int i = 0; i < static_cast< int >( Profile::Count ); ++i )
		SetParamElementInfo( PT_PROFILE, i, kProfileNames[ i ], static_cast< float >( i ) );
	SetParamInfo( PT_PELLET, "Pellet", FF_TYPE_EVENT, false );
	SetParamInfof( PT_FEED, "Feed", FF_TYPE_STANDARD );
	SetParamInfof( PT_CLIP_HEATS, "Clip Heats", FF_TYPE_STANDARD );
	SetParamInfof( PT_FUEL, "Fuel", FF_TYPE_STANDARD );

	SetParamInfof( PT_FIELD, "Field", FF_TYPE_STANDARD );
	SetParamInfof( PT_GUIDE_FIELD, "Guide Field", FF_TYPE_STANDARD );
	SetOptionParamInfo( PT_POLES, "Poles", kPoleOptions, params[ PT_POLES ] );
	for( int i = 0; i < kPoleOptions; ++i )
		SetParamElementInfo( PT_POLES, i, kPoleNames[ i ], static_cast< float >( i ) );
	SetParamInfof( PT_COIL_RADIUS, "Coil Radius", FF_TYPE_STANDARD );
	SetParamInfof( PT_COIL_SPIN, "Coil Spin", FF_TYPE_STANDARD );
	SetParamInfof( PT_CURVATURE, "Curvature", FF_TYPE_STANDARD );
	SetParamInfo( PT_QUENCH, "Quench", FF_TYPE_BOOLEAN, false );
	SetOptionParamInfo( PT_BOUNDARY, "Boundary", static_cast< int >( Boundary::Count ), params[ PT_BOUNDARY ] );
	for( int i = 0; i < static_cast< int >( Boundary::Count ); ++i )
		SetParamElementInfo( PT_BOUNDARY, i, kBoundaryNames[ i ], static_cast< float >( i ) );

	SetParamInfof( PT_SPEED, "Speed", FF_TYPE_STANDARD );
	SetParamInfof( PT_RESISTIVITY, "Resistivity", FF_TYPE_STANDARD );
	SetParamInfof( PT_COOLING, "Cooling", FF_TYPE_STANDARD );
	SetOptionParamInfo( PT_DETAIL, "Detail", kDetailCount, params[ PT_DETAIL ] );
	for( int i = 0; i < kDetailCount; ++i )
		SetParamElementInfo( PT_DETAIL, i, kDetailNames[ i ], static_cast< float >( i ) );
	SetParamInfof( PT_DRIVE, "Drive", FF_TYPE_STANDARD );
	SetParamInfof( PT_DRIVE_SCALE, "Drive Scale", FF_TYPE_STANDARD );

	//An FFT buffer: Resolume shows it as an audio-source picker and writes the
	//spectrum into it every frame.
	SetBufferParamInfo( PT_AUDIO, "Audio", audio::kBins, FF_USAGE_FFT );
	for( int i = 0; i < audio::kBins; ++i )
		SetParamElementInfo( PT_AUDIO, i, "", 0.0f );
	SetParamInfof( PT_AUDIO_HEAT, "Audio Heat", FF_TYPE_STANDARD );
	SetParamInfof( PT_AUDIO_PELLETS, "Audio Pellets", FF_TYPE_STANDARD );

	SetParamInfof( PT_EXPOSURE, "Exposure", FF_TYPE_STANDARD );
	SetParamInfof( PT_TINT, "Temperature Tint", FF_TYPE_STANDARD );
	SetOptionParamInfo( PT_RAMP, "Ramp", static_cast< int >( Ramp::Count ), params[ PT_RAMP ] );
	for( int i = 0; i < static_cast< int >( Ramp::Count ); ++i )
		SetParamElementInfo( PT_RAMP, i, kRampNames[ i ], static_cast< float >( i ) );
	SetParamInfof( PT_GLOW, "Glow", FF_TYPE_STANDARD );
	SetParamInfof( PT_FIELD_LINES, "Field Lines", FF_TYPE_STANDARD );
	SetParamInfof( PT_LINE_COUNT, "Field Line Count", FF_TYPE_STANDARD );
	SetOptionParamInfo( PT_VIEW, "View", static_cast< int >( View::Count ), params[ PT_VIEW ] );
	for( int i = 0; i < static_cast< int >( View::Count ); ++i )
		SetParamElementInfo( PT_VIEW, i, kViewNames[ i ], static_cast< float >( i ) );
	SetParamInfof( PT_MIX, "Mix", FF_TYPE_STANDARD );

	//Groups: SetParamGroup collapses runs of consecutive same-group ids.
	SetParamGroup( PT_PRESET, "Preset" );
	for( unsigned int id = PT_IGNITE; id <= PT_FUEL; ++id )
		SetParamGroup( id, "Ball" );
	for( unsigned int id = PT_FIELD; id <= PT_BOUNDARY; ++id )
		SetParamGroup( id, "Bottle" );
	for( unsigned int id = PT_SPEED; id <= PT_DRIVE_SCALE; ++id )
		SetParamGroup( id, "Plasma" );
	for( unsigned int id = PT_AUDIO; id <= PT_AUDIO_PELLETS; ++id )
		SetParamGroup( id, "Audio" );
	for( unsigned int id = PT_EXPOSURE; id <= PT_MIX; ++id )
		SetParamGroup( id, "Light" );

	SetParamInfo( PT_ABOUT_TEXT, "About", FF_TYPE_TEXT, stoatworks::about::defaultText() );
	{
		FFUInt32 aboutId = PT_ABOUT_TEXT + 1;
		for( const auto& b : stoatworks::about::buttons() )
			SetParamInfo( aboutId++, b.label, FF_TYPE_EVENT, false );
	}
	for( unsigned int id = PT_ABOUT_TEXT; id < PT_COUNT; ++id )
		SetParamGroup( id, "About" );
}

//---------------------------------------------------------------------------
FFResult ContainmentPlugin::InitGL( const FFGLViewportStruct* vp )
{
	diag::init();
	diag::info( std::string( "GL vendor=" ) + glStringOrUnknown( GL_VENDOR ) + " renderer="
	            + glStringOrUnknown( GL_RENDERER ) + " version=" + glStringOrUnknown( GL_VERSION ) );

	for( int i = 0; i < static_cast< int >( Program::Count ); ++i )
	{
		ProgramSource source = SourceFor( static_cast< Program >( i ) );
		if( !g_mutationFind.empty() )
		{
			const size_t at = source.fragment.find( g_mutationFind );
			if( at != std::string::npos )
				source.fragment.replace( at, g_mutationFind.size(), g_mutationReplace );
		}
		if( programs[ i ].Compile( source.vertex, source.fragment ) )
			continue;

		//Returning FF_FAIL is invisible to the operator: the plugin simply does
		//nothing. These lines are the only record of which pass it was.
		diag::error( std::string( "the " ) + source.name + " shader failed to compile - the plugin will do nothing" );
		FFGLLog::LogToHost( "Containment: shader failed to compile" );
		DeInitGL();
		return FF_FAIL;
	}

	if( !quad.Initialise() )
	{
		diag::error( "quad geometry failed to initialise" );
		DeInitGL();
		return FF_FAIL;
	}

	diag::info( "initialised" );
	return CFFGLPlugin::InitGL( vp );
}

//---------------------------------------------------------------------------
ContainmentPlugin::Model ContainmentPlugin::CurrentModel() const
{
	Model m;
	m.field     = FieldFromParam( P( PT_FIELD ) );
	m.guide     = GuideFieldFromParam( P( PT_GUIDE_FIELD ) ) * m.field;
	m.poles     = kPoleCounts[ optionIndex( P( PT_POLES ), kPoleOptions ) ];
	m.spin      = CoilSpinFromParam( P( PT_COIL_SPIN ) );
	m.curvature = CurvatureFromParam( P( PT_CURVATURE ) );
	m.boundary  = test.boundary >= 0 ? test.boundary
	                                 : optionIndex( P( PT_BOUNDARY ), static_cast< int >( Boundary::Count ) );
	m.eta        = ResistivityFromParam( P( PT_RESISTIVITY ) );
	m.cooling    = CoolingFromParam( P( PT_COOLING ) );
	m.drive      = DriveFromParam( P( PT_DRIVE ) );
	m.driveScale = DriveScaleFromParam( P( PT_DRIVE_SCALE ) );

	//The coils never sit inside the frame.
	const double halfDiagonal = 0.5 * std::sqrt( grid.lx * grid.lx + grid.ly * grid.ly );
	m.coilRadius = std::max( static_cast< double >( CoilRadiusFromParam( P( PT_COIL_RADIUS ) ) ), 1.05 * halfDiagonal );

	//Ball X / Y are fractions of the FRAME, which starts at the margin.
	m.ballX      = ( grid.ox + std::clamp( P( PT_BALL_X ), 0.0f, 1.0f ) * grid.fx ) * grid.dx;
	m.ballY      = ( grid.oy + std::clamp( P( PT_BALL_Y ), 0.0f, 1.0f ) * grid.fy ) * grid.dx;
	m.ballRadius = BallSizeFromParam( P( PT_BALL_SIZE ) );
	m.pressure   = 0.5 * TemperatureFromParam( P( PT_TEMPERATURE ) );
	return m;
}

//---------------------------------------------------------------------------
bool ContainmentPlugin::EnsureBuffers( int width, int height, const Grid& wanted )
{
	const bool reset = wanted.nx != grid.nx || wanted.ny != grid.ny;
	if( reset )
		ignitePending = true;

	bool ok = state[ 0 ].Ensure( wanted.nx, wanted.ny ) && state[ 1 ].Ensure( wanted.nx, wanted.ny )
	          && star.Ensure( wanted.nx, wanted.ny ) && fluxX.Ensure( wanted.nx + 1, wanted.ny )
	          && fluxY.Ensure( wanted.nx, wanted.ny + 1 );

	//The reduction: kReduceBlock-square blocks, until one more level would
	//leave at most a few hundred texels for the clock pass to finish.
	if( reset || reduction.empty() )
	{
		for( PassBuffer& b : reduction )
			b.Destroy();
		reduction.clear();
		int w = wanted.nx, h = wanted.ny;
		do
		{
			w = ( w + kReduceBlock - 1 ) / kReduceBlock;
			h = ( h + kReduceBlock - 1 ) / kReduceBlock;
			reduction.emplace_back();
		} while( w * h > 16 );
	}
	{
		int w = wanted.nx, h = wanted.ny;
		for( PassBuffer& b : reduction )
		{
			w  = ( w + kReduceBlock - 1 ) / kReduceBlock;
			h  = ( h + kReduceBlock - 1 ) / kReduceBlock;
			ok = ok && b.Ensure( w, h, GL_RGBA32F, PassBuffer::Sampling::Nearest );
		}
	}
	for( PassBuffer& b : clock )
		ok = ok && b.Ensure( 2, 1, GL_RGBA32F, PassBuffer::Sampling::Nearest );

	//The light, on the grid; the glow on a copy reduced to <= kGlowCells.
	//The light is the frame's, not the margin's.
	ok = ok && emission.Ensure( wanted.fx, wanted.fy, GL_RGBA32F, PassBuffer::Sampling::Linear );
	int gw = wanted.fx, gh = wanted.fy, reductions = 0;
	while( std::min( gw, gh ) > kGlowCells && gw % 2 == 0 && gh % 2 == 0 )
	{
		gw /= 2;
		gh /= 2;
		++reductions;
	}
	if( static_cast< int >( glowReduce.size() ) != reductions )
	{
		for( PassBuffer& b : glowReduce )
			b.Destroy();
		glowReduce.clear();
		glowReduce.resize( static_cast< size_t >( reductions ) );
	}
	{
		int w = wanted.fx, h = wanted.fy;
		for( PassBuffer& b : glowReduce )
		{
			w /= 2;
			h /= 2;
			ok = ok && b.Ensure( w, h, GL_RGBA32F, PassBuffer::Sampling::Nearest );
		}
	}
	ok = ok && glowTemp.Ensure( gw, gh, GL_RGBA32F, PassBuffer::Sampling::Nearest );
	for( PassBuffer& b : glow )
		ok = ok && b.Ensure( gw, gh, GL_RGBA32F, PassBuffer::Sampling::Linear );
	glowWidth  = gw;
	glowHeight = gh;

	//The multigrid for A_z: halve while both sides stay even and >= 8.
	if( reset || levels.empty() )
	{
		for( Level& level : levels )
		{
			level.solution[ 0 ].Destroy();
			level.solution[ 1 ].Destroy();
			level.rhs.Destroy();
			level.residual.Destroy();
		}
		levels.clear();
		int w = wanted.nx, h = wanted.ny;
		double spacing = wanted.dx;
		while( true )
		{
			Level level;
			level.nx      = w;
			level.ny      = h;
			level.spacing = spacing;
			levels.push_back( std::move( level ) );
			if( w % 2 != 0 || h % 2 != 0 || std::min( w, h ) < 16 )
				break;
			w /= 2;
			h /= 2;
			spacing *= 2.0;
		}
		havePotential = false;
	}
	for( Level& level : levels )
	{
		ok = ok && level.solution[ 0 ].Ensure( level.nx, level.ny, GL_RGBA32F, PassBuffer::Sampling::Linear );
		ok = ok && level.solution[ 1 ].Ensure( level.nx, level.ny, GL_RGBA32F, PassBuffer::Sampling::Linear );
		ok = ok && level.rhs.Ensure( level.nx, level.ny, GL_RGBA32F, PassBuffer::Sampling::Nearest );
		ok = ok && level.residual.Ensure( level.nx, level.ny, GL_RGBA32F, PassBuffer::Sampling::Nearest );
	}

	if( !ok )
		return false;
	grid = wanted;
	return true;
}

//---------------------------------------------------------------------------
void ContainmentPlugin::SetStateUniforms( GLuint p ) const
{
	uniform2i( p, "GridSize", grid.nx, grid.ny );
	uniform1f( p, "Dx", static_cast< float >( grid.dx ) );
	uniform2i( p, "FrameOrigin", grid.ox, grid.oy );
	uniform2i( p, "FrameCells", grid.fx, grid.fy );
	uniform1f( p, "Gamma", test.gamma );
	const Model m = CurrentModel();
	uniform1i( p, "BoundaryMode", m.boundary );
	uniform1i( p, "UseFloors", test.floors ? 1 : 0 );
	uniform1f( p, "RhoFloor", kDensityFloor );
	uniform1f( p, "PFloor", kPressureFloor );
	uniform1f( p, "AmbientDensity", test.backgroundDensity );
	uniform1f( p, "AmbientPressure", test.backgroundPressure );
	uniform1i( p, "Solver", test.solver );
	uniform1i( p, "UseGLM", test.glm ? 1 : 0 );

	float data[ kMaxCoils * 3 ] = {};
	for( int k = 0; k < coils.count; ++k )
	{
		data[ k * 3 + 0 ] = static_cast< float >( coils.x[ k ] );
		data[ k * 3 + 1 ] = static_cast< float >( coils.y[ k ] );
		data[ k * 3 + 2 ] = static_cast< float >( coils.current[ k ] );
	}
	uniform1i( p, "CoilCount", coils.count );
	glUniform3fv( location( p, "CoilData" ), kMaxCoils, data );
	uniform1f( p, "GuideBz", static_cast< float >( coils.guide ) );

	uniform1f( p, "Curvature", static_cast< float >( m.curvature ) );
	//The reference temperature: the ball's as ignited, p0 / rho0.
	uniform1f( p, "CurvatureTemp", test.uniformGravity ? 0.0f : static_cast< float >( m.pressure ) );
	if( test.gravityCentre )
		uniform2f( p, "GravityCentre", static_cast< float >( test.gravityX ), static_cast< float >( test.gravityY ) );
	else
		uniform2f( p, "GravityCentre", static_cast< float >( m.ballX ), static_cast< float >( m.ballY ) );
	uniform1f( p, "GravityCore",
	           static_cast< float >( test.gravityCore >= 0.0 ? test.gravityCore : 0.5 * m.ballRadius ) );

	float modes[ kDriveModes * 4 ] = {};
	int count                      = 0;
	if( m.drive > 0.0 )
		for( const Drive::Mode& mode : drive.Modes() )
		{
			modes[ count * 4 + 0 ] = mode.kx;
			modes[ count * 4 + 1 ] = mode.ky;
			modes[ count * 4 + 2 ] = mode.a;
			modes[ count * 4 + 3 ] = mode.b;
			++count;
		}
	uniform1i( p, "DriveCount", count );
	uniform1f( p, "DriveWindow", static_cast< float >( kDriveWindow * m.ballRadius ) );
	glUniform4fv( location( p, "DriveModes" ), kDriveModes, modes );
}

void ContainmentPlugin::DrawInto( const StateBuffer& target )
{
	glBindFramebuffer( GL_FRAMEBUFFER, target.Framebuffer() );
	glViewport( 0, 0, target.Width(), target.Height() );
	quad.Draw();
}

//---------------------------------------------------------------------------
void ContainmentPlugin::Ignite( GLuint input, float maxU, float maxV )
{
	const Model m   = CurrentModel();
	const GLuint p  = Id( Program::Ignite );
	glUseProgram( p );
	SetStateUniforms( p );
	samplers( p, { "InputTexture" } );
	uniform2f( p, "MaxUV", maxU, maxV );
	uniform2f( p, "BallCentre", static_cast< float >( m.ballX ), static_cast< float >( m.ballY ) );
	uniform1f( p, "BallRadius", static_cast< float >( m.ballRadius ) );
	uniform1i( p, "ProfileKind", optionIndex( P( PT_PROFILE ), static_cast< int >( Profile::Count ) ) );
	uniform1f( p, "BallPressure", static_cast< float >( m.pressure ) );
	uniform1f( p, "ClipHeats", ClipHeatsFromParam( P( PT_CLIP_HEATS ) ) );
	uniform1f( p, "BgDensity", test.backgroundDensity );
	uniform1f( p, "BgPressure", test.backgroundPressure );
	{
		BoundTextures bound( { input } );
		DrawInto( state[ current ] );
	}
	glUseProgram( 0 );

	drive.Reset();
	speedKnown    = false;
	havePotential = false;
}

void ContainmentPlugin::Sources( GLuint input, float maxU, float maxV, double hostDt, double simDt, int pellets,
                                 bool always )
{
	const Model m  = CurrentModel();
	const float feed = FeedFromParam( P( PT_FEED ) );
	//Feed is the fraction per frame at 60 fps; frame-rate independent.
	const float fraction = 1.0f - static_cast< float >( std::pow( 1.0 - std::min( feed, 0.999f ), hostDt * 60.0 ) );
	//Heating, in p0 per tau_A: the clip's light through Feed, and the audio.
	const float clipHeat  = static_cast< float >( 2.0 * ClipHeatsFromParam( P( PT_CLIP_HEATS ) ) * feed
                                                 * m.pressure * simDt );
	const float audioHeat = static_cast< float >( AudioHeatFromParam( P( PT_AUDIO_HEAT ) ) * analyser.Level()
	                                              * m.pressure * simDt );
	const float fuel = static_cast< float >( 1.0 - std::exp( -FuelFromParam( P( PT_FUEL ) ) * simDt ) );
	if( fraction <= 0.0f && clipHeat <= 0.0f && audioHeat <= 0.0f && pellets == 0 && fuel <= 0.0f && !always
	    && !carryCoils )
		return;

	const GLuint p = Id( Program::Sources );
	glUseProgram( p );
	SetStateUniforms( p );
	samplers( p, { "StateA", "StateB", "StateC", "StateD", "InputTexture" } );
	uniform2f( p, "MaxUV", maxU, maxV );
	uniform1f( p, "FeedFraction", fraction );
	uniform1f( p, "ClipHeatRate", clipHeat );
	uniform1f( p, "AudioHeatRate", audioHeat );
	uniform1i( p, "PelletCount", pellets );
	uniform1f( p, "FuelFraction", always ? 0.0f : fuel );
	uniform2f( p, "BallCentre", static_cast< float >( m.ballX ), static_cast< float >( m.ballY ) );
	uniform1f( p, "BallRadius", static_cast< float >( m.ballRadius ) );
	uniform1i( p, "ProfileKind", optionIndex( P( PT_PROFILE ), static_cast< int >( Profile::Count ) ) );
	uniform1f( p, "BallPressure", static_cast< float >( m.pressure ) );
	uniform1f( p, "BgDensity", test.backgroundDensity );
	uniform1f( p, "BgPressure", test.backgroundPressure );
	uniform2f( p, "PelletCentre", static_cast< float >( m.ballX ), static_cast< float >( m.ballY ) );
	uniform1f( p, "PelletDensity", kPelletDensity );
	uniform1f( p, "PelletRadius", kPelletRadius );
	float data[ kMaxCoils * 3 ] = {};
	for( int k = 0; k < previousCoils.count; ++k )
	{
		data[ k * 3 + 0 ] = static_cast< float >( previousCoils.x[ k ] );
		data[ k * 3 + 1 ] = static_cast< float >( previousCoils.y[ k ] );
		data[ k * 3 + 2 ] = static_cast< float >( previousCoils.current[ k ] );
	}
	uniform1i( p, "CarryCoils", carryCoils && !always ? 1 : 0 );
	uniform1i( p, "PrevCoilCount", previousCoils.count );
	glUniform3fv( location( p, "PrevCoilData" ), kMaxCoils, data );
	uniform1f( p, "PrevGuideBz", static_cast< float >( previousCoils.guide ) );
	{
		const StateBuffer& s = state[ current ];
		BoundTextures bound( { s.Texture( 0 ), s.Texture( 1 ), s.Texture( 2 ), s.Texture( 3 ), input } );
		DrawInto( state[ 1 - current ] );
	}
	glUseProgram( 0 );
	current = 1 - current;
}

//---------------------------------------------------------------------------
void ContainmentPlugin::Clock( int substep, int substeps, double target )
{
	//The reduction.
	{
		const GLuint p = Id( Program::Reduce );
		glUseProgram( p );
		samplers( p, { "Source" } );
		GLuint source = state[ current ].Texture( 2 );
		int sw = grid.nx, sh = grid.ny;
		bool fromState = true;
		uniform1i( p, "Block", kReduceBlock );
		for( PassBuffer& level : reduction )
		{
			uniform2i( p, "SourceSize", sw, sh );
			uniform1i( p, "FromState", fromState ? 1 : 0 );
			BoundTextures bound( { source } );
			glBindFramebuffer( GL_FRAMEBUFFER, level.GetGLID() );
			glViewport( 0, 0, level.GetWidth(), level.GetHeight() );
			quad.Draw();
			source    = level.TextureID();
			sw        = level.GetWidth();
			sh        = level.GetHeight();
			fromState = false;
		}
	}

	//The clock.
	const Model m     = CurrentModel();
	const GLuint p    = Id( Program::Clock );
	const int next    = 1 - clockIndex;
	glUseProgram( p );
	samplers( p, { "Reduced", "PrevClock" } );
	uniform2i( p, "ReducedSize", reduction.back().GetWidth(), reduction.back().GetHeight() );
	uniform1i( p, "Substep", substep );
	uniform1i( p, "Substeps", substeps );
	uniform1f( p, "Target", static_cast< float >( target ) );
	uniform1f( p, "Courant", kCourant );
	uniform1f( p, "Dx", static_cast< float >( grid.dx ) );
	uniform1f( p, "Eta", static_cast< float >( m.eta ) );
	{
		BoundTextures bound( { reduction.back().TextureID(), clock[ clockIndex ].TextureID() } );
		glBindFramebuffer( GL_FRAMEBUFFER, clock[ next ].GetGLID() );
		glViewport( 0, 0, 2, 1 );
		quad.Draw();
	}
	glUseProgram( 0 );
	clockIndex = next;
}

void ContainmentPlugin::Substep()
{
	const Model m        = CurrentModel();
	const StateBuffer& s = state[ current ];

	//MUSCL-Hancock's half step.
	{
		const GLuint p = Id( Program::Predict );
		glUseProgram( p );
		SetStateUniforms( p );
		samplers( p, { "StateA", "StateB", "StateC", "StateD", "ClockTexture" } );
		BoundTextures bound(
			{ s.Texture( 0 ), s.Texture( 1 ), s.Texture( 2 ), s.Texture( 3 ), clock[ clockIndex ].TextureID() } );
		DrawInto( star );
	}

	//Every face's flux, once.
	{
		const GLuint p = Id( Program::Flux );
		glUseProgram( p );
		SetStateUniforms( p );
		samplers( p, { "StateA", "StateB", "StateC", "StateD", "StarA", "StarB", "StarC", "StarD", "ClockTexture" } );
		uniform1f( p, "Eta", static_cast< float >( m.eta ) );
		BoundTextures bound( { s.Texture( 0 ), s.Texture( 1 ), s.Texture( 2 ), s.Texture( 3 ), star.Texture( 0 ),
		                       star.Texture( 1 ), star.Texture( 2 ), star.Texture( 3 ),
		                       clock[ clockIndex ].TextureID() } );
		uniform1i( p, "Direction", 0 );
		DrawInto( fluxX );
		uniform1i( p, "Direction", 1 );
		DrawInto( fluxY );
	}

	//The update.
	{
		const GLuint p = Id( Program::Update );
		glUseProgram( p );
		SetStateUniforms( p );
		samplers( p, { "StateA", "StateB", "StateC", "StateD", "StarA", "StarB", "ClockTexture", "FluxXA", "FluxXB",
		               "FluxXC", "FluxXD", "FluxYA", "FluxYB", "FluxYC", "FluxYD" } );
		uniform1f( p, "CoolingRate", static_cast< float >( m.cooling ) );
		uniform1f( p, "CoolingFloorT", test.backgroundPressure / test.backgroundDensity );
		uniform1f( p, "GLMAlpha", kGLMAlpha );
		uniform1i( p, "UseEntropy", test.entropy ? 1 : 0 );
		uniform1f( p, "EntropySwitch", kEntropySwitch );
		const float efolds = test.spongeEFolds >= 0.0f ? test.spongeEFolds : kSpongeEFolds;
		uniform1f( p, "SpongeEFolds", grid.ox > 0 ? efolds : 0.0f );
		uniform1f( p, "SpongePower", kSpongePower );
		BoundTextures bound( { s.Texture( 0 ), s.Texture( 1 ), s.Texture( 2 ), s.Texture( 3 ), star.Texture( 0 ),
		                       star.Texture( 1 ), clock[ clockIndex ].TextureID(), fluxX.Texture( 0 ), fluxX.Texture( 1 ),
		                       fluxX.Texture( 2 ), fluxX.Texture( 3 ), fluxY.Texture( 0 ), fluxY.Texture( 1 ),
		                       fluxY.Texture( 2 ), fluxY.Texture( 3 ) } );
		DrawInto( state[ 1 - current ] );
	}
	glUseProgram( 0 );
	current = 1 - current;
	++substepsTaken;
}

void ContainmentPlugin::ReadClock( float out[ 4 ] )
{
	float both[ 8 ] = {};
	glBindTexture( GL_TEXTURE_2D, clock[ clockIndex ].TextureID() );
	glGetTexImage( GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, both );
	glBindTexture( GL_TEXTURE_2D, 0 );
	for( int i = 0; i < 4; ++i )
		out[ i ] = both[ i ];
	lastEntropyCells = both[ 4 ];
}

void ContainmentPlugin::Advance( double target, bool synchronous )
{
	if( target <= 0.0 )
		return;

	//The plan needs the fastest signal speed. After an Ignite (or with no
	//history) it is measured now, synchronously; otherwise last frame's is
	//used, with 10% to spare -- and if the plasma sped up since, each substep
	//still takes only the CFL step, so the frame covers less time.
	if( synchronous || !speedKnown )
	{
		Clock( 0, 1, 0.0 );
		float values[ 4 ] = {};
		ReadClock( values );
		lastSpeed  = 2.0f * values[ 2 ];
		speedKnown = true;
	}
	const double stable = kCourant * grid.dx / std::max( static_cast< double >( lastSpeed ), 1e-12 );
	double step         = stable;
	const Model m       = CurrentModel();
	if( m.eta > 0.0 )
		step = std::min( step, 0.2 * grid.dx * grid.dx / m.eta );
	const int wanted    = static_cast< int >( std::ceil( 1.1 * target / step ) );
	frameSubsteps       = std::clamp( wanted, 1, kMaxSubsteps );

	for( int k = 0; k < frameSubsteps; ++k )
	{
		Clock( k, frameSubsteps, target );
		Substep();
	}
	lastTarget   = target;
	clockPending = true;
}

int ContainmentPlugin::StepForTest( double duration, int limit )
{
	int taken       = 0;
	double remaining = duration;
	//The stirring moves on once per call, as it does once per frame.
	{
		const Model m = CurrentModel();
		drive.Advance( duration, m.driveScale, m.drive );
	}
	while( remaining > duration * 1e-7 && taken < limit )
	{
		//Plan exactly, from the speeds as they stand.
		Clock( 0, 1, 0.0 );
		float values[ 4 ] = {};
		ReadClock( values );
		if( !std::isfinite( values[ 2 ] ) || values[ 2 ] > 1e30f )
			return -1;
		lastSpeed  = 2.0f * values[ 2 ];
		speedKnown = true;
		const Model m = CurrentModel();
		double step   = kCourant * grid.dx / std::max( static_cast< double >( lastSpeed ), 1e-12 );
		if( m.eta > 0.0 )
			step = std::min( step, 0.2 * grid.dx * grid.dx / m.eta );
		const int n = std::clamp( static_cast< int >( std::ceil( 1.05 * remaining / step ) ), 1, kMaxSubsteps );
		for( int k = 0; k < n; ++k )
		{
			Clock( k, n, remaining );
			Substep();
		}
		taken += n;
		ReadClock( values );
		floorHits += values[ 3 ];
		entropyCells += lastEntropyCells;
		if( !std::isfinite( values[ 1 ] ) || values[ 1 ] <= 0.0f )
			return -1;
		remaining -= values[ 1 ];
		simTime += values[ 1 ];
	}
	//The last update's flags are only counted by a reduction after it.
	{
		Clock( 0, 1, 0.0 );
		float values[ 4 ] = {};
		ReadClock( values );
		floorHits += values[ 3 ];
	}
	return remaining > duration * 1e-6 ? -1 : taken;
}

void ContainmentPlugin::LoadStateForTest( const std::vector< float >& a, const std::vector< float >& b,
                                          const std::vector< float >& c, const std::vector< float >& d )
{
	const size_t n = static_cast< size_t >( grid.nx ) * grid.ny * 4;
	const std::vector< float >* data[ 4 ] = { &a, &b, &c, &d };
	for( int i = 0; i < 4; ++i )
	{
		if( data[ i ]->size() != n )
			return;
		glBindTexture( GL_TEXTURE_2D, state[ current ].Texture( i ) );
		glTexSubImage2D( GL_TEXTURE_2D, 0, 0, 0, grid.nx, grid.ny, GL_RGBA, GL_FLOAT, data[ i ]->data() );
	}
	glBindTexture( GL_TEXTURE_2D, 0 );

	//The loaded state carries no signal speeds; a sources pass that adds
	//nothing writes them. (Its clip reads are of the state's own texture,
	//harmlessly: nothing is fed.)
	const float feed = P( PT_FEED );
	params[ PT_FEED ] = 0.0f;
	Sources( state[ current ].Texture( 0 ), 1.0f, 1.0f, 0.0, 0.0, 0, true );
	params[ PT_FEED ] = feed;

	ignitePending = false;
	speedKnown    = false;
	havePotential = false;
	drive.Reset();
}

//---------------------------------------------------------------------------
void ContainmentPlugin::Emission()
{
	const GLuint p = Id( Program::Emission );
	glUseProgram( p );
	SetStateUniforms( p );
	samplers( p, { "StateA", "StateB", "StateC", "StateD" } );
	uniform1i( p, "View", optionIndex( P( PT_VIEW ), static_cast< int >( View::Count ) ) );
	uniform1f( p, "Tint", std::clamp( P( PT_TINT ), 0.0f, 1.0f ) );
	uniform1i( p, "RampKind", optionIndex( P( PT_RAMP ), static_cast< int >( Ramp::Count ) ) );
	uniform1f( p, "TempReference", static_cast< float >( CurrentModel().pressure ) );
	uniform1f( p, "FieldReference", static_cast< float >( CurrentModel().field ) );
	const StateBuffer& s = state[ current ];
	{
		BoundTextures bound( { s.Texture( 0 ), s.Texture( 1 ), s.Texture( 2 ), s.Texture( 3 ) } );
		glBindFramebuffer( GL_FRAMEBUFFER, emission.GetGLID() );
		glViewport( 0, 0, grid.fx, grid.fy );
		quad.Draw();
	}
	glUseProgram( 0 );
}

void ContainmentPlugin::Glow()
{
	glowFraction = GlowFromParam( P( PT_GLOW ) );
	if( glowFraction <= 0.0f )
		return;

	//Reduce, exactly: every 2x2 block to its mean.
	GLuint source = emission.TextureID();
	{
		const GLuint p = Id( Program::Downsample );
		glUseProgram( p );
		samplers( p, { "Source" } );
		for( PassBuffer& level : glowReduce )
		{
			BoundTextures bound( { source } );
			glBindFramebuffer( GL_FRAMEBUFFER, level.GetGLID() );
			glViewport( 0, 0, level.GetWidth(), level.GetHeight() );
			quad.Draw();
			source = level.TextureID();
		}
	}

	//Three Gaussians, cascaded: each stage blurs the last by the difference in
	//variance, so stage i is the emission under a Gaussian of kGlowSigma[ i ].
	const GLuint p = Id( Program::Blur );
	glUseProgram( p );
	samplers( p, { "Source" } );
	uniform2i( p, "Size", glowWidth, glowHeight );
	double previous = 0.0;
	std::vector< float > weights;
	for( int stage = 0; stage < kGlowStages; ++stage )
	{
		//The glow's copy spans the frame, whose height is 1.
		const double sigma = kGlowSigma[ stage ] * glowHeight;
		const double step  = std::sqrt( std::max( sigma * sigma - previous * previous, 1e-6 ) );
		previous           = sigma;
		const int taps     = gaussianWeights( step, weights );
		uniform1i( p, "Taps", taps );
		glUniform1fv( location( p, "Weights" ), kMaxBlurTaps + 1, weights.data() );

		struct Pass
		{
			GLuint from;
			GLuint into;
			int dx, dy;
		};
		const Pass passes[ 2 ] = { { source, glowTemp.GetGLID(), 1, 0 },
			                       { glowTemp.TextureID(), glow[ stage ].GetGLID(), 0, 1 } };
		for( const Pass& pass : passes )
		{
			uniform2i( p, "Direction", pass.dx, pass.dy );
			BoundTextures bound( { pass.from } );
			glBindFramebuffer( GL_FRAMEBUFFER, pass.into );
			glViewport( 0, 0, glowWidth, glowHeight );
			quad.Draw();
		}
		source = glow[ stage ].TextureID();
	}
	glUseProgram( 0 );
}

void ContainmentPlugin::Potential()
{
	const GLuint residual = Id( Program::Residual );
	const GLuint smooth   = Id( Program::Smooth );
	const GLuint prolong  = Id( Program::Prolong );

	float data[ kMaxCoils * 3 ] = {};
	for( int k = 0; k < coils.count; ++k )
	{
		data[ k * 3 + 0 ] = static_cast< float >( coils.x[ k ] );
		data[ k * 3 + 1 ] = static_cast< float >( coils.y[ k ] );
		data[ k * 3 + 2 ] = static_cast< float >( coils.current[ k ] );
	}
	auto common = [ & ]( GLuint p, const Level& level, bool finest ) {
		uniform2i( p, "Size", level.nx, level.ny );
		uniform1f( p, "Spacing", static_cast< float >( level.spacing ) );
		uniform1i( p, "VacuumGhosts", finest ? 1 : 0 );
		uniform1i( p, "CoilCount", coils.count );
		glUniform3fv( location( p, "CoilData" ), kMaxCoils, data );
		samplers( p, { "Solution", "Rhs", "StateB", "Coarse" } );
	};
	auto target = [ & ]( PassBuffer& buffer ) {
		glBindFramebuffer( GL_FRAMEBUFFER, buffer.GetGLID() );
		glViewport( 0, 0, buffer.GetWidth(), buffer.GetHeight() );
	};
	auto relax = [ & ]( Level& level, bool finest, int sweeps ) {
		glUseProgram( smooth );
		common( smooth, level, finest );
		for( int i = 0; i < sweeps; ++i )
		{
			BoundTextures bound( { level.solution[ level.current ].TextureID(), level.rhs.TextureID() } );
			target( level.solution[ 1 - level.current ] );
			quad.Draw();
			level.current = 1 - level.current;
		}
	};

	Level& top = levels[ 0 ];
	//The right-hand side, -J_z, from the field as it stands.
	glUseProgram( residual );
	common( residual, top, true );
	uniform1i( residual, "Mode", 0 );
	{
		//Every sampler the program declares gets a real texture, even the ones
		//this mode never reads: a unit left at 0 under a declared sampler is
		//logged by the driver as an unloadable texture.
		BoundTextures bound( { top.residual.TextureID(), top.residual.TextureID(), state[ current ].Texture( 1 ),
		                       top.residual.TextureID() } );
		target( top.rhs );
		quad.Draw();
	}
	//A cold start: the coils' own potential is the best guess there is.
	if( !havePotential )
	{
		top.solution[ 0 ].Clear();
		top.solution[ 1 ].Clear();
		relax( top, true, 20 );
		havePotential = true;
	}

	//One V-cycle a frame, seeded from the last frame's answer.
	const int count = static_cast< int >( levels.size() );
	for( int l = 0; l < count; ++l )
	{
		Level& level = levels[ l ];
		if( l > 0 )
		{
			level.solution[ 0 ].Clear();
			level.solution[ 1 ].Clear();
			level.current = 0;
		}
		const bool coarsest = l == count - 1;
		relax( level, l == 0, coarsest ? 40 : 2 );
		if( coarsest )
			break;

		glUseProgram( residual );
		common( residual, level, l == 0 );
		uniform1i( residual, "Mode", 1 );
		{
			BoundTextures bound(
				{ level.solution[ level.current ].TextureID(), level.rhs.TextureID(), state[ current ].Texture( 1 ) } );
			target( level.residual );
			quad.Draw();
		}
		const GLuint down = Id( Program::Downsample );
		glUseProgram( down );
		samplers( down, { "Source" } );
		{
			BoundTextures bound( { level.residual.TextureID() } );
			target( levels[ l + 1 ].rhs );
			quad.Draw();
		}
	}
	for( int l = count - 2; l >= 0; --l )
	{
		Level& level = levels[ l ];
		glUseProgram( prolong );
		common( prolong, level, l == 0 );
		{
			const Level& coarse = levels[ l + 1 ];
			BoundTextures bound( { level.solution[ level.current ].TextureID(), level.rhs.TextureID(),
			                       state[ current ].Texture( 1 ), coarse.solution[ coarse.current ].TextureID() } );
			target( level.solution[ 1 - level.current ] );
			quad.Draw();
		}
		level.current = 1 - level.current;
		relax( level, l == 0, 2 );
	}
	glUseProgram( 0 );
}

//---------------------------------------------------------------------------
void ContainmentPlugin::UpdateClock()
{
	const double wallNow = wallSeconds();
	if( wallStart < 0.0 )
		wallStart = wallNow;
	const double raw = hostTime;

	//Resolume has been seen sending both seconds and milliseconds. Vote on it
	//against the wall clock for a few frames, then stop asking (rosette).
	if( clockScale == 0.0 && raw >= 0.0 && lastRawTime >= 0.0 && lastWallTime >= 0.0 )
	{
		const double hostDelta = raw - lastRawTime;
		const double wallDelta = wallNow - lastWallTime;
		if( hostDelta > 0.0 && wallDelta >= 0.0005 )
		{
			const double ratio = hostDelta / wallDelta;
			if( ratio > 0.1 && ratio < 10.0 )
				++secondsVotes;
			else if( ratio > 100.0 && ratio < 10000.0 )
				++millisVotes;
			if( secondsVotes >= kClockVotes || millisVotes >= kClockVotes )
			{
				clockScale  = millisVotes > secondsVotes ? 0.001 : 1.0;
				settledJump = true;
			}
		}
	}
	if( raw >= 0.0 )
		lastRawTime = raw;
	lastWallTime = wallNow;

	//Resolume's clock overflows a float (the fleet measured ~5e8 ms), so only
	//DIFFERENCES of it are ever used, in double.
	now = ( raw >= 0.0 && clockScale != 0.0 ) ? raw * clockScale : wallNow - wallStart;
	if( ++clockFrames == 60 )
		diag::info( "host clock at frame 60: raw=" + std::to_string( raw ) + " scale=" + std::to_string( clockScale ) );
}

//---------------------------------------------------------------------------
FFResult ContainmentPlugin::ProcessOpenGL( ProcessOpenGLStruct* pgl )
{
	if( pgl == nullptr || pgl->numInputTextures < 1 || pgl->inputTextures[ 0 ] == nullptr )
		return FF_FAIL;

	const FFGLTextureStruct& source = *( pgl->inputTextures[ 0 ] );
	const int width                 = static_cast< int >( source.Width );
	const int height                = static_cast< int >( source.Height );
	if( width <= 0 || height <= 0 )
		return FF_FAIL;

	ScopedGLState restore;
	const GLint* hostViewport = restore.saved.viewport;
	glDisable( GL_BLEND );

	//-------------------------------------------------------------------
	// Time. Frame-relative and clamped.
	//-------------------------------------------------------------------
	UpdateClock();
	if( settledJump )
	{
		lastNow     = -1.0;
		settledJump = false;
	}
	const double hostDt = lastNow >= 0.0 ? std::clamp( now - lastNow, 0.0, kMaxFrameDelta ) : 0.0;
	lastNow             = now;

	float bins[ audio::kBins ] = {};
	int binCount               = 0;
	if( const ParamInfo* info = FindParamInfo( PT_AUDIO ) )
	{
		binCount = static_cast< int >( std::min< size_t >( info->elements.size(), audio::kBins ) );
		for( int i = 0; i < binCount; ++i )
			bins[ i ] = info->elements[ static_cast< size_t >( i ) ].value;
	}
	audio::Settings settings;
	settings.sensitivity = std::clamp( P( PT_AUDIO_PELLETS ), 0.0f, 1.0f );
	analyser.Update( bins, binCount, static_cast< float >( hostDt ), settings );

	//Last frame's clock, read now rather than then, so the GPU has had a
	//frame to finish it.
	if( clockPending )
	{
		float values[ 4 ] = {};
		ReadClock( values );
		lastSpeed = 2.0f * values[ 2 ];
		floorHits += values[ 3 ];
		entropyCells += lastEntropyCells;
		clockPending = false;
		if( values[ 1 ] < 0.999f * static_cast< float >( lastTarget ) )
		{
			++cappedFrames;
			if( ++framesSinceCapLog > 600 )
			{
				framesSinceCapLog = 0;
				diag::warn( "the substep cap (" + std::to_string( kMaxSubsteps ) + ") bit: a frame covered "
				            + std::to_string( values[ 1 ] ) + " of " + std::to_string( lastTarget )
				            + " tau_A -- simulated time is running slow (" + std::to_string( cappedFrames )
				            + " frames so far)" );
			}
		}
		else
			++framesSinceCapLog;
	}

	//-------------------------------------------------------------------
	// The grid, and every allocation, before anything binds a texture.
	//-------------------------------------------------------------------
	const int detail  = optionIndex( P( PT_DETAIL ), kDetailCount );
	//Open is a window onto a bigger bottle: a margin of simulated cells all
	//round the frame, never shown, holding the absorbing layer.
	const int boundaryNow = test.boundary >= 0 ? test.boundary
	                                           : optionIndex( P( PT_BOUNDARY ), static_cast< int >( Boundary::Count ) );
	const float marginFraction = test.marginFraction >= 0.0f ? test.marginFraction : kMarginFraction;
	const int margin  = boundaryNow == static_cast< int >( Boundary::Open ) && !test.legacyOpen
	                        ? static_cast< int >( std::lround( marginFraction * kDetailCells[ detail ] ) )
	                        : 0;
	const Grid wanted = gridRasterW > 0 && gridRasterH > 0
	                        ? ChooseGrid( gridRasterW, gridRasterH, kDetailCells[ detail ], margin )
	                        : ChooseGrid( width, height, kDetailCells[ detail ], margin );
	if( !EnsureBuffers( width, height, wanted ) )
	{
		diag::error( "could not allocate the plasma's buffers" );
		return FF_FAIL;
	}

	//-------------------------------------------------------------------
	// The coils, as of this frame.
	//-------------------------------------------------------------------
	const Model m      = CurrentModel();
	const double simDt = hostDt * SpeedFromParam( P( PT_SPEED ) );//tau_A = 1
	spinAngle += m.spin * simDt;
	spinAngle = std::remainder( spinAngle, 2.0 * kPi );
	//The quench: the coils' current decays with kQuenchTime, and comes back
	//the same way when Quench is let go.
	const double wantStrength = P( PT_QUENCH ) >= 0.5f ? 0.0 : 1.0;
	coilStrength += ( wantStrength - coilStrength ) * ( 1.0 - std::exp( -simDt / kQuenchTime ) );
	previousCoils = coils;
	coils = MakeCoils( m.poles, m.coilRadius, spinAngle, m.field, m.guide, coilStrength, 0.5 * grid.lx, 0.5 * grid.ly );
	//With a conducting vessel the coils' change is carried through it (the
	//sources pass). Any change: a turn, a quench, a new Field or Coil Radius
	//-- but not a new pole count, whose coils are not the same coils.
	//Carried through the volume for Open as well: imposing the turn only at
	//the edges left the interior lagging and piled the mismatch up as current
	//sheets where the edges face the coils, which flung the plasma (|v| 18 at
	//4.5 tau_A on the default look).
	carryCoils = ( m.boundary == static_cast< int >( Boundary::Wall ) || m.boundary == static_cast< int >( Boundary::Open ) )
	             && previousCoils.count == coils.count && !ignitePending && grid.nx > 0;
	drive.Advance( simDt, m.driveScale, m.drive );

	const FFGLTexCoords maxCoords = GetMaxGLTexCoords( source );

	//-------------------------------------------------------------------
	// Events, then the plasma.
	//-------------------------------------------------------------------
	bool synchronous = false;
	if( ignitePending )
	{
		Ignite( source.Handle, maxCoords.s, maxCoords.t );
		ignitePending = false;
		synchronous   = true;
	}
	int pellets   = pelletPresses;
	pelletPresses = 0;
	if( P( PT_AUDIO_PELLETS ) > 0.0f && analyser.Fired() )
		++pellets;
	Sources( source.Handle, maxCoords.s, maxCoords.t, hostDt, simDt, pellets );
	if( pellets > 0 )
		synchronous = true;

	Advance( simDt, synchronous );
	simTime += simDt;

	//-------------------------------------------------------------------
	// The light.
	//-------------------------------------------------------------------
	Emission();
	const int view = optionIndex( P( PT_VIEW ), static_cast< int >( View::Count ) );
	if( view == 0 )
		Glow();
	const float lines   = std::clamp( P( PT_FIELD_LINES ), 0.0f, 1.0f );
	float lineSpacing   = 0.0f;
	if( view == 0 && lines > 0.0f && coils.count > 0 )
	{
		Potential();
		//Count contours across the bottle's flux, centre to rim: the largest
		//|A - A_centre| round the circle of radius 0.5. (Not in a gap: there
		//B is radial, so dA/dr = -B_phi = 0 and a multipole's A is 0 -- the
		//first version measured the span there, got ~0, and drew the lines
		//infinitely close.)
		const double cx = 0.5 * grid.lx, cy = 0.5 * grid.ly;
		const double centre = VacuumPotential( coils, cx, cy );
		double span = 0.0;
		for( int k = 0; k < 360; ++k )
		{
			const double a = 2.0 * kPi * k / 360.0;
			span = std::max( span, std::abs( VacuumPotential( coils, cx + 0.5 * std::cos( a ), cy + 0.5 * std::sin( a ) ) - centre ) );
		}
		lineSpacing = static_cast< float >( span / LineCountFromParam( P( PT_LINE_COUNT ) ) );
	}

	//-------------------------------------------------------------------
	// Composite, into the host's framebuffer and viewport.
	//-------------------------------------------------------------------
	glBindFramebuffer( GL_FRAMEBUFFER, pgl->HostFBO );
	glViewport( hostViewport[ 0 ], hostViewport[ 1 ], hostViewport[ 2 ], hostViewport[ 3 ] );
	{
		const GLuint p = Id( Program::Composite );
		glUseProgram( p );
		samplers( p, { "InputTexture", "Emission", "Glow0", "Glow1", "Glow2", "Potential" } );
		uniform2f( p, "MaxUV", maxCoords.s, maxCoords.t );
		const bool exact = source.HardwareWidth == source.Width && source.HardwareHeight == source.Height
		                   && hostViewport[ 2 ] == width && hostViewport[ 3 ] == height && hostViewport[ 0 ] == 0
		                   && hostViewport[ 1 ] == 0;
		uniform1i( p, "ExactInput", exact ? 1 : 0 );
		uniform1i( p, "View", view );
		const double ev = ExposureFromParam( P( PT_EXPOSURE ) );
		uniform1f( p, "Gain", static_cast< float >( std::pow( 2.0, ev ) / kEmissionReference ) );
		const float glowAmount = view == 0 ? glowFraction : 0.0f;
		uniform1f( p, "GlowAmount", glowAmount );
		//The light moved into the glare leaves the core. The wrong model, for
		//--negative, leaves the core whole and ADDS the glare: a bloom.
		uniform1f( p, "CoreWeight", test.additiveGlow ? 1.0f : 1.0f - glowAmount );
		uniform3f( p, "GlowShare", kGlowShare[ 0 ], kGlowShare[ 1 ], kGlowShare[ 2 ] );
		uniform1f( p, "LineAmount", lines );
		uniform1f( p, "LineSpacing", lineSpacing );
		//The potential covers the margin; the frame's uv maps into it.
		glUniform4f( location( p, "PotentialMap" ), static_cast< float >( grid.fx ) / grid.nx,
		             static_cast< float >( grid.fy ) / grid.ny, static_cast< float >( grid.ox ) / grid.nx,
		             static_cast< float >( grid.oy ) / grid.ny );
		uniform1f( p, "MixAmount", std::clamp( P( PT_MIX ), 0.0f, 1.0f ) );
		const GLuint potential = havePotential ? PotentialTextureID() : emission.TextureID();
		BoundTextures bound( { source.Handle, emission.TextureID(), glow[ 0 ].TextureID(), glow[ 1 ].TextureID(),
		                       glow[ 2 ].TextureID(), potential } );
		quad.Draw();
		glUseProgram( 0 );
	}

	return FF_SUCCESS;
}

//---------------------------------------------------------------------------
FFResult ContainmentPlugin::DeInitGL()
{
	for( ffglex::FFGLShader& program : programs )
		program.FreeGLResources();
	quad.Release();

	for( StateBuffer& s : state )
		s.Destroy();
	star.Destroy();
	fluxX.Destroy();
	fluxY.Destroy();
	for( PassBuffer& b : reduction )
		b.Destroy();
	reduction.clear();
	for( PassBuffer& b : clock )
		b.Destroy();
	emission.Destroy();
	for( PassBuffer& b : glowReduce )
		b.Destroy();
	glowReduce.clear();
	glowTemp.Destroy();
	for( PassBuffer& b : glow )
		b.Destroy();
	for( Level& level : levels )
	{
		level.solution[ 0 ].Destroy();
		level.solution[ 1 ].Destroy();
		level.rhs.Destroy();
		level.residual.Destroy();
	}
	levels.clear();

	grid          = Grid {};
	ignitePending = true;
	speedKnown    = false;
	clockPending  = false;
	havePotential = false;
	return FF_SUCCESS;
}

//---------------------------------------------------------------------------
FFResult ContainmentPlugin::SetTime( double time )
{
	hostTime = time;
	return FF_SUCCESS;
}

char* ContainmentPlugin::GetTextParameter( unsigned int index )
{
	if( index == PT_ABOUT_TEXT )
	{
		static const std::string text = stoatworks::about::textParam( 0 );
		return const_cast< char* >( text.c_str() );
	}
	return CFFGLPlugin::GetTextParameter( index );
}

FFResult ContainmentPlugin::SetTextParameter( unsigned int index, const char* value )
{
	if( index == PT_ABOUT_TEXT )
		return FF_SUCCESS;
	return CFFGLPlugin::SetTextParameter( index, value );
}

FFResult ContainmentPlugin::SetFloatParameter( unsigned int index, float value )
{
	if( index >= PT_COUNT )
		return FF_FAIL;
	if( index >= PT_ABOUT_TEXT )
		return stoatworks::about::handleParam( index - PT_ABOUT_TEXT, value ) ? FF_SUCCESS : FF_FAIL;

	//Events act on the rising edge: a host sends 1 then 0 and may restate
	//either.
	const bool down = value >= 0.5f;
	if( index == PT_IGNITE )
	{
		if( down && !igniteHeld )
			ignitePending = true;
		igniteHeld = down;
	}
	else if( index == PT_PRESET )
	{
		//A new bottle wants a fresh ball.
		if( std::lround( value ) != std::lround( params[ PT_PRESET ] ) )
			ignitePending = true;
	}
	else if( index == PT_PELLET )
	{
		if( down && !pelletHeld )
			++pelletPresses;
		pelletHeld = down;
	}

	params[ index ] = value;
	return FF_SUCCESS;
}

float ContainmentPlugin::GetFloatParameter( unsigned int index )
{
	if( index >= PT_COUNT )
		return 0.0f;
	return params[ index ];
}

} // namespace containment
