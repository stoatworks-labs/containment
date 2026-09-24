/**
    cttest -- render Containment offline, and measure what its plasma is doing.

    It drives the REAL plugin class, through the same ProcessOpenGL a host
    calls, on a synthetic 60 fps clock. A test that exercises a
    reimplementation tests the reimplementation.

        cttest --out /tmp/frame.png       the card, as a plasma, after --frames
        cttest --list                     every parameter and its default
        cttest --preset FILE              apply 'Name=value' lines, like --set
        cttest --pipe                     raw frames in, raw frames out
        cttest --film N                   N frames of the card, raw frames out
        cttest --dump-shaders DIR         every program exactly as compiled

    `--script` is the fleet's cue format, `frame  Parameter Name  value`, held
    before the first key and interpolated between -- so a button press is
    three keys (29 Ignite 0 / 30 Ignite 1 / 31 Ignite 0).

    The claims, one flag each (see README / AGENTS.md for every tolerance):

        --briowu     the Brio-Wu shock tube, along x and along y, at two Details
        --alfven     a circularly polarised Alfven wave: speed, order, rotation
        --conserve   Wall: mass and energy to float round-off over 600 frames
        --divb       div B stays small on a turbulent run; GLM off must fail
        --balance    the diamagnetic bubble: p + B^2/2 flat, Bz inside, edge, beta = 1
        --rt         magnetic Rayleigh-Taylor growth rates, three cases
        --cusp       leaks at the N cusp angles, N = 4, 6, 8, and following the spin
        --frozen     Bz / rho is carried with the fluid
        --quench     the free expansion's front and the coils' decay
        --resist     field diffuses into a static plasma on L^2 / eta
        --floors     the floors fire rarely; without them the plasma blows up
        --still      Mix 0 is the identity, bit for bit
        --glow       the glare moves light and makes none
        --state      the host's GL state comes back as it went in
        --mutation   one character changed in a shipped shader must fail a check
        --negative   every check above against a deliberately wrong model
        --bench      ms/frame at 720p, 1080p and 4K for each Detail
*/

#include "Checks.h"
#include "Harness.h"

#include <algorithm>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <unistd.h>

using namespace cttest;

namespace
{
struct NamedParameter
{
	std::string name;
	unsigned int index;
	float value;
	std::string kind;
};

const char* kindName( unsigned int type )
{
	switch( type )
	{
	case FF_TYPE_BOOLEAN: return "bool";
	case FF_TYPE_EVENT: return "event";
	case FF_TYPE_XPOS: return "xpos";
	case FF_TYPE_YPOS: return "ypos";
	case FF_TYPE_OPTION: return "option";
	case FF_TYPE_BUFFER: return "buffer";
	case FF_TYPE_STANDARD: return "standard";
	case FF_TYPE_TEXT: return "text";
	default: return "other";
	}
}

std::vector< NamedParameter > listParameters( ContainmentPlugin& plugin )
{
	std::vector< NamedParameter > list;
	for( unsigned int i = 0; i < PT_COUNT; ++i )
	{
		const char* const name = plugin.GetParamName( i );
		list.push_back( { name ? name : "?", i, plugin.GetFloatParameter( i ), kindName( plugin.GetParamType( i ) ) } );
	}
	return list;
}

bool applySetting( ContainmentPlugin& plugin, const std::string& assignment, std::string& error )
{
	const size_t equals = assignment.find( '=' );
	if( equals == std::string::npos )
	{
		error = "expected Name=Value";
		return false;
	}
	const std::string name  = assignment.substr( 0, equals );
	const std::string value = assignment.substr( equals + 1 );
	for( const NamedParameter& parameter : listParameters( plugin ) )
		if( parameter.name == name )
		{
			plugin.SetFloatParameter( parameter.index, std::strtof( value.c_str(), nullptr ) );
			return true;
		}
	error = "no parameter called '" + name + "'";
	return false;
}

/// A preset: one `Name=value` per line, `#` comments.
bool loadPreset( const std::string& path, std::vector< std::string >& settings )
{
	std::ifstream file( path );
	if( !file )
		return false;
	std::string line;
	while( std::getline( file, line ) )
	{
		const size_t hash = line.find( '#' );
		if( hash != std::string::npos )
			line.erase( hash );
		while( !line.empty() && ( line.back() == ' ' || line.back() == '\t' || line.back() == '\r' ) )
			line.pop_back();
		size_t start = line.find_first_not_of( " \t" );
		if( start == std::string::npos )
			continue;
		settings.push_back( line.substr( start ) );
	}
	return true;
}

using Track = std::vector< std::pair< int, float > >;

std::map< std::string, Track > loadScript( const std::string& path, std::string& error )
{
	std::map< std::string, Track > tracks;
	std::ifstream file( path );
	if( !file )
	{
		error = "cannot open " + path;
		return tracks;
	}
	std::string line;
	int lineNumber = 0;
	while( std::getline( file, line ) )
	{
		++lineNumber;
		const size_t hash = line.find( '#' );
		if( hash != std::string::npos )
			line.erase( hash );
		std::istringstream in( line );
		int frame = 0;
		if( !( in >> frame ) )
			continue;
		std::vector< std::string > words;
		std::string word;
		while( in >> word )
			words.push_back( word );
		if( words.size() < 2 )
		{
			error = path + ":" + std::to_string( lineNumber ) + ": expected `frame Parameter Name value`";
			return {};
		}
		const float value = std::strtof( words.back().c_str(), nullptr );
		words.pop_back();
		std::string name = words.front();
		for( size_t i = 1; i < words.size(); ++i )
			name += " " + words[ i ];
		tracks[ name ].emplace_back( frame, value );
	}
	for( auto& entry : tracks )
		std::sort( entry.second.begin(), entry.second.end() );
	return tracks;
}

float valueAt( const Track& track, int frame )
{
	if( track.empty() )
		return 0.0f;
	if( frame <= track.front().first )
		return track.front().second;
	if( frame >= track.back().first )
		return track.back().second;
	for( size_t i = 0; i + 1 < track.size(); ++i )
	{
		const auto& a = track[ i ];
		const auto& b = track[ i + 1 ];
		if( frame >= a.first && frame <= b.first )
		{
			if( b.first == a.first )
				return b.second;
			const float t = static_cast< float >( frame - a.first ) / static_cast< float >( b.first - a.first );
			return a.second + ( b.second - a.second ) * t;
		}
	}
	return track.back().second;
}

int runPipe( int width, int height, const std::string& scriptPath, int filmFrames, bool beat,
             const std::vector< std::string >& settings )
{
	Rig rig;
	if( !rig.Init( width, height ) )
		return 1;
	if( beat )
		rig.feed = AudioFeed::Pulses;
	for( const std::string& setting : settings )
	{
		std::string error;
		if( !applySetting( rig.plugin, setting, error ) )
		{
			std::fprintf( stderr, "--set %s: %s\n", setting.c_str(), error.c_str() );
			return 2;
		}
	}

	std::map< unsigned int, Track > automation;
	if( !scriptPath.empty() )
	{
		std::string error;
		const std::map< std::string, Track > tracks = loadScript( scriptPath, error );
		if( !error.empty() )
		{
			std::fprintf( stderr, "%s\n", error.c_str() );
			return 2;
		}
		const std::vector< NamedParameter > known = listParameters( rig.plugin );
		for( const auto& entry : tracks )
		{
			bool found = false;
			for( const NamedParameter& parameter : known )
				if( parameter.name == entry.first )
				{
					automation[ parameter.index ] = entry.second;
					found                         = true;
				}
			if( !found )
			{
				std::fprintf( stderr, "script names '%s', which is not a parameter (try --list)\n", entry.first.c_str() );
				return 2;
			}
		}
	}

	std::vector< unsigned char > in( static_cast< size_t >( width ) * height * 4 );
	Floats picture( in.size() );
	for( int index = 0; filmFrames < 0 || index < filmFrames; ++index )
	{
		if( filmFrames < 0 )
		{
			size_t filled = 0;
			while( filled < in.size() )
			{
				const ssize_t got = read( STDIN_FILENO, in.data() + filled, in.size() - filled );
				if( got <= 0 )
					break;
				filled += static_cast< size_t >( got );
			}
			if( filled < in.size() )
				break;
			for( int y = 0; y < height; ++y )
				for( int x = 0; x < width * 4; ++x )
					picture[ static_cast< size_t >( height - 1 - y ) * width * 4 + x ] =
						in[ static_cast< size_t >( y ) * width * 4 + x ] / 255.0f;
			rig.Upload( picture );
		}
		for( const auto& track : automation )
			rig.plugin.SetFloatParameter( track.first, valueAt( track.second, index ) );
		if( !rig.Render( 1 ) )
			return 1;

		const Floats out = rig.Output();
		std::vector< unsigned char > bytes( in.size() );
		for( int y = 0; y < height; ++y )
			for( int x = 0; x < width * 4; ++x )
				bytes[ static_cast< size_t >( y ) * width * 4 + x ] = static_cast< unsigned char >( std::lround(
					std::clamp( out[ static_cast< size_t >( height - 1 - y ) * width * 4 + x ], 0.0f, 1.0f ) * 255.0f ) );
		size_t written = 0;
		while( written < bytes.size() )
		{
			const ssize_t put = write( STDOUT_FILENO, bytes.data() + written, bytes.size() - written );
			if( put <= 0 )
			{
				std::fprintf( stderr, "cttest --pipe: stdout closed after %d frames (%s)\n", index,
				              put < 0 ? std::strerror( errno ) : "no progress" );
				return 1;
			}
			written += static_cast< size_t >( put );
		}
	}
	return 0;
}

int dumpShaders( const std::string& directory )
{
	int n = 0;
	for( const ProgramSource& source : AllSources() )
	{
		std::ofstream v( directory + "/" + source.name + ".vert" );
		v << source.vertex;
		std::ofstream f( directory + "/" + source.name + ".frag" );
		f << source.fragment;
		//A directory that is not there must not read as "wrote 14 programs".
		if( !v || !f )
		{
			std::fprintf( stderr, "cannot write %s into %s\n", source.name, directory.c_str() );
			return 1;
		}
		++n;
	}
	std::printf( "%d\n", n );
	return 0;
}
} // namespace

//---------------------------------------------------------------------------
int main( int argc, char** argv )
{
	std::string outPath = "/tmp/containment.png";
	std::vector< std::string > settings;
	int width = 1280, height = 720, frames = 180;
	std::vector< int > igniteFrames, pelletFrames;
	bool beat = false;
	std::string mode, scriptPath, dumpDirectory;
	int filmFrames = -1;
	bool sizeGiven = false, allowNoGL = false;

	for( int i = 1; i < argc; ++i )
	{
		const std::string argument = argv[ i ];
		const bool hasNext         = i + 1 < argc;
		if( argument == "--help" || argument == "-h" )
		{
			std::printf(
				"cttest -- render Containment offline and measure its plasma\n\n"
				"  --out PATH        render the card and write it here\n"
				"  --size WxH        render size (default 1280x720)\n"
				"  --frames N        frames at 60 fps before reading back (default 180)\n"
				"  --ignite N        press Ignite on frame N. Repeatable.\n"
				"  --pellet N        press Pellet on frame N. Repeatable.\n"
				"  --beat            feed a beat every half second into the Audio buffer\n"
				"  --set \"Name=V\"    set a parameter by its display name. Repeatable.\n"
				"  --preset FILE     apply a file of Name=V lines\n"
				"  --list            print every parameter and its default, then exit\n"
				"  --pipe            raw RGBA frames on stdin, raw RGBA frames on stdout\n"
				"  --film N          N frames of the card, raw RGBA frames on stdout\n"
				"  --script PATH     parameter cues for --pipe/--film: 'frame Name value'\n"
				"  --dump-shaders D  write every program, as compiled, into directory D\n\n"
				"  --offline         every check that needs no GL context, and their negative controls\n"
				"  --allow-no-gl     with a GL check: SKIP loudly, not FAIL, when no context can be made\n"
				"  --size WxH        with a check: every rig renders at this raster instead of its own\n\n"
				"  --briowu --alfven --conserve --divb --balance --rt --cusp --frozen --quench\n"
				"  --resist --floors --still --glow --state --open --presets --mutation --negative\n"
				"  --reference --vacuum --names (offline) --bench\n" );
			return 0;
		}
		else if( argument == "--out" && hasNext )
			outPath = argv[ ++i ];
		else if( argument == "--set" && hasNext )
			settings.push_back( argv[ ++i ] );
		else if( argument == "--preset" && hasNext )
		{
			const std::string path = argv[ ++i ];
			if( !loadPreset( path, settings ) )
			{
				std::fprintf( stderr, "cannot read preset %s\n", path.c_str() );
				return 2;
			}
		}
		else if( argument == "--frames" && hasNext )
			frames = std::atoi( argv[ ++i ] );
		else if( argument == "--ignite" && hasNext )
			igniteFrames.push_back( std::atoi( argv[ ++i ] ) );
		else if( argument == "--pellet" && hasNext )
			pelletFrames.push_back( std::atoi( argv[ ++i ] ) );
		else if( argument == "--beat" )
			beat = true;
		else if( argument == "--pipe" )
			mode = "pipe";
		else if( argument == "--film" && hasNext )
		{
			mode       = "pipe";
			filmFrames = std::max( 1, std::atoi( argv[ ++i ] ) );
		}
		else if( argument == "--script" && hasNext )
			scriptPath = argv[ ++i ];
		else if( argument == "--dump-shaders" && hasNext )
		{
			mode          = "dump";
			dumpDirectory = argv[ ++i ];
		}
		else if( argument == "--list" )
			mode = "list";
		else if( argument == "--size" && hasNext )
		{
			const std::string value = argv[ ++i ];
			const size_t cross      = value.find( 'x' );
			if( cross != std::string::npos )
			{
				width  = std::atoi( value.substr( 0, cross ).c_str() );
				height = std::atoi( value.substr( cross + 1 ).c_str() );
			}
			if( cross == std::string::npos || width < 8 || height < 8 || width > 8192 || height > 8192 )
			{
				std::fprintf( stderr, "--size wants WxH, each 8..8192 (got '%s')\n", value.c_str() );
				return 2;
			}
			sizeGiven = true;
		}
		else if( argument == "--allow-no-gl" )
			allowNoGL = true;
		else if( argument.rfind( "--", 0 ) == 0 )
			mode = argument.substr( 2 );
		else
		{
			std::fprintf( stderr, "unknown argument '%s' (try --help)\n", argument.c_str() );
			return 2;
		}
	}

	if( mode == "list" )
	{
		ContainmentPlugin plugin;
		std::printf( "%-3s %-18s %-9s %s\n", "id", "name", "kind", "default" );
		for( const NamedParameter& parameter : listParameters( plugin ) )
			std::printf( "%-3u %-18s %-9s %.4f\n", parameter.index, parameter.name.c_str(), parameter.kind.c_str(),
			             parameter.value );
		return 0;
	}
	if( mode == "dump" )
		return dumpShaders( dumpDirectory );

	//--offline is every check that needs no GL context, listed HERE, in the
	//one place that knows which those are: a GitHub macOS runner cannot make
	//an accelerated context, and a workflow that listed them itself would go
	//stale the first time one was added.
	{
		const Perturb none;
		const std::pair< const char*, int ( * )( const Perturb& ) > offline[] = {
			{ "names", RunNames }, { "presets", RunPresets }, { "reference", RunReference }, { "vacuum", RunVacuum } };
		for( const auto& check : offline )
			if( mode == check.first )
			{
				check.second( none );
				return Verdict();
			}
		if( mode == "offline" )
		{
			for( const auto& check : offline )
				check.second( none );
			const int checks = Verdict();
			const int negatives = RunNegative( true );
			std::printf( "\n  OFFLINE: nothing above drew a pixel through a GL driver. The solver, the bottle, the\n"
			             "  light, every GL check and every GL negative control were NOT run -- tools/verify.sh\n"
			             "  runs them on a GPU. In CI the shaders were only compiled, by glslc.\n" );
			return checks != 0 || negatives != 0 ? 1 : 0;
		}
	}

	//A check's rigs render at --size when it is given (verify.sh's 320x180
	//pass); the pipe and the plain render always use it.
	if( sizeGiven && mode != "pipe" && !mode.empty() && mode != "stats" )
	{
		g_rasterW = width;
		g_rasterH = height;
	}

	CGLContextObj context = CreateContext();
	if( context == nullptr )
	{
		if( allowNoGL )
		{
			std::printf( "  SKIP  could not create an OpenGL 4.1 core context, accelerated or software. --%s was NOT "
			             "run.\n",
			             mode.empty() ? "out" : mode.c_str() );
			return 0;
		}
		std::fprintf( stderr, "could not create an OpenGL 4.1 core context\n" );
		return 1;
	}
	if( !mode.empty() && mode != "pipe" && mode != "stats" && mode != "bench" )
		std::printf( "GL %s / %s%s\n", glGetString( GL_VERSION ), glGetString( GL_RENDERER ),
		             g_rasterW > 0 ? fmt( ", every rig at %dx%d", g_rasterW, g_rasterH ).c_str() : ", each check at its own raster" );

	const Perturb none;
	const std::pair< const char*, int ( * )( const Perturb& ) > checks[] = {
		{ "briowu", RunBrioWu }, { "alfven", RunAlfven }, { "conserve", RunConserve }, { "divb", RunDivB },
		{ "balance", RunBalance }, { "rt", RunRT },       { "cusp", RunCusp },         { "frozen", RunFrozen },
		{ "quench", RunQuench }, { "resist", RunResist }, { "floors", RunFloors },     { "still", RunStill },
		{ "glow", RunGlow },     { "state", RunState },   { "mutation", RunMutation }, { "equilibrium", RunEquilibrium },
		{ "open", RunOpen },
	};

	int result = -1;
	for( const auto& check : checks )
		if( mode == check.first )
		{
			check.second( none );
			result = Verdict();
		}

	if( result < 0 )
	{
		if( mode == "pipe" )
		{
			//A reader that hangs up (`| head -c 1`, ffmpeg dying) must end the
			//take with exit 1 and a message, not SIGPIPE's silent 141: write()
			//then fails and runPipe says so.
			std::signal( SIGPIPE, SIG_IGN );
			result = runPipe( width, height, scriptPath, filmFrames, beat, settings );
		}
		else if( mode == "negative" )
			result = RunNegative();
		else if( mode == "bench" )
			result = RunBench( [ & ]( Rig& rig ) {
				for( const std::string& setting : settings )
				{
					std::string error;
					if( !applySetting( rig.plugin, setting, error ) )
						std::fprintf( stderr, "--set %s: %s\n", setting.c_str(), error.c_str() );
				}
			} );
		else if( mode == "stats" )
		{
			//A development aid: the state's ranges after --frames, with --set.
			Rig rig;
			result = rig.Init( width, height ) ? 0 : 1;
			for( const std::string& setting : settings )
			{
				std::string error;
				applySetting( rig.plugin, setting, error );
			}
			for( int f = 0; f < frames && result == 0; f += std::max( 1, frames / 6 ) )
			{
				rig.Render( std::max( 1, frames / 6 ) );
				const Snapshot s = Snapshot::Take( rig.plugin );
				double rmin = 1e30, rmax = -1e30, pmin = 1e30, pmax = -1e30, bmax = 0, vmax = 0;
				int bad = 0;
				for( int j = 0; j < s.ny; ++j )
					for( int i = 0; i < s.nx; ++i )
					{
						const double r = s.Rho( i, j ), p = s.P( i, j );
						if( !std::isfinite( r ) || !std::isfinite( p ) )
						{
							++bad;
							continue;
						}
						rmin = std::min( rmin, r ); rmax = std::max( rmax, r );
						pmin = std::min( pmin, p ); pmax = std::max( pmax, p );
						bmax = std::max( bmax, std::sqrt( s.Bx( i, j ) * s.Bx( i, j ) + s.By( i, j ) * s.By( i, j ) + s.Bz( i, j ) * s.Bz( i, j ) ) );
						vmax = std::max( vmax, std::sqrt( s.U( i, j ) * s.U( i, j ) + s.V( i, j ) * s.V( i, j ) ) );
					}
				const Floats e = ReadTexture( rig.plugin.EmissionTextureID(), rig.plugin.EmissionWidth(), rig.plugin.EmissionHeight() );
				double emax = 0;
				for( size_t k = 3; k < e.size(); k += 4 )
					emax = std::max( emax, static_cast< double >( e[ k ] ) );
				if( const char* dump = std::getenv( "CT_FLAGS" ) )
				{
					Floats img( static_cast< size_t >( s.nx ) * s.ny * 4 );
					for( int j = 0; j < s.ny; ++j )
						for( int i = 0; i < s.nx; ++i )
						{
							const size_t o = s.At( i, j );
							img[ o + 0 ] = static_cast< float >( std::min( 1.0, std::hypot( s.U( i, j ), s.V( i, j ) ) / 3.0 ) );
							img[ o + 1 ] = static_cast< float >( std::clamp( ( std::log10( std::max( s.P( i, j ), 1e-9 ) ) + 6.0 ) / 6.0, 0.0, 1.0 ) );
							img[ o + 2 ] = static_cast< float >( std::clamp( s.Rho( i, j ), 0.0, 1.0 ) );
							img[ o + 3 ] = 1.0f;
						}
					WritePng( std::string( dump ) + std::to_string( rig.frame ) + ".png", s.nx, s.ny, img );
				}
				std::printf( "frame %4d t=%.4f grid %dx%d rho [%.4g %.4g] p [%.4g %.4g] |B|max %.3g |v|max %.3g emax %.3g nonfinite %d substeps %lld floors %.0f capped %d\n",
				             rig.frame, rig.plugin.SimTime(), s.nx, s.ny, rmin, rmax, pmin, pmax, bmax, vmax, emax, bad,
				             rig.plugin.SubstepsTaken(), rig.plugin.FloorHits(), rig.plugin.CappedFrames() );
			}
		}
		else if( !mode.empty() )
		{
			std::fprintf( stderr, "unknown mode --%s (try --help)\n", mode.c_str() );
			result = 2;
		}
		else
		{
			Rig rig;
			result = rig.Init( width, height ) ? 0 : 1;
			for( const std::string& setting : settings )
			{
				std::string error;
				if( result == 0 && !applySetting( rig.plugin, setting, error ) )
				{
					std::fprintf( stderr, "--set %s: %s\n", setting.c_str(), error.c_str() );
					result = 2;
				}
			}
			if( beat )
				rig.feed = AudioFeed::Pulses;
			for( int f = 0; f < std::max( frames, 1 ) && result == 0; ++f )
			{
				if( std::find( igniteFrames.begin(), igniteFrames.end(), f ) != igniteFrames.end() )
					rig.Press( PT_IGNITE );
				if( std::find( pelletFrames.begin(), pelletFrames.end(), f ) != pelletFrames.end() )
					rig.Press( PT_PELLET );
				if( !rig.Render( 1 ) )
					result = 1;
			}
			if( result == 0 )
			{
				if( WritePng( outPath, width, height, rig.Output() ) )
					std::printf( "wrote %s -- %dx%d, %d frames (%.3f tau_A), %lld substeps, %d capped frames\n",
					             outPath.c_str(), width, height, frames, rig.plugin.SimTime(),
					             rig.plugin.SubstepsTaken(), rig.plugin.CappedFrames() );
				else
					result = 1;
			}
		}
	}

	CGLSetCurrentContext( nullptr );
	CGLDestroyContext( context );
	return result;
}
