#include "Harness.h"

#include <cstdlib>
#include <string>

#include <zlib.h>

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstring>

namespace cttest
{
int g_failures = 0;
int g_rasterW   = 0;
int g_rasterH   = 0;

void ChooseRaster( int& width, int& height )
{
	if( g_rasterW > 0 && g_rasterH > 0 )
	{
		width  = g_rasterW;
		height = g_rasterH;
	}
}

std::string fmt( const char* format, ... )
{
	char buffer[ 2048 ];
	va_list args;
	va_start( args, format );
	std::vsnprintf( buffer, sizeof( buffer ), format, args );
	va_end( args );
	return buffer;
}

void Check( bool condition, const std::string& message )
{
	std::printf( "  %s  %s\n", condition ? "ok  " : "FAIL", message.c_str() );
	std::fflush( stdout );
	if( !condition )
		++g_failures;
}

void Say( const char* format, ... )
{
	va_list args;
	va_start( args, format );
	std::vprintf( format, args );
	va_end( args );
	std::fflush( stdout );
}

int Verdict()
{
	std::printf( "\n  %s\n", g_failures == 0 ? "PASS" : "FAIL" );
	return g_failures == 0 ? 0 : 1;
}

//---------------------------------------------------------------------------
CGLContextObj CreateContext()
{
	const CGLPixelFormatAttribute accelerated[] = {
		kCGLPFAOpenGLProfile, static_cast< CGLPixelFormatAttribute >( kCGLOGLPVersion_GL4_Core ),
		kCGLPFAAccelerated,
		kCGLPFAColorSize, static_cast< CGLPixelFormatAttribute >( 24 ),
		kCGLPFAAlphaSize, static_cast< CGLPixelFormatAttribute >( 8 ),
		static_cast< CGLPixelFormatAttribute >( 0 )
	};
	const CGLPixelFormatAttribute software[] = {
		kCGLPFAOpenGLProfile, static_cast< CGLPixelFormatAttribute >( kCGLOGLPVersion_GL4_Core ),
		kCGLPFAColorSize, static_cast< CGLPixelFormatAttribute >( 24 ),
		kCGLPFAAlphaSize, static_cast< CGLPixelFormatAttribute >( 8 ),
		static_cast< CGLPixelFormatAttribute >( 0 )
	};
	//CT_RENDERER=software: Apple's generic float renderer, a different
	//rasteriser with GL 4.1's legal float slack -- what a GPU-less CI runner
	//would give -- on this Mac, so "would this hold on another rasteriser?"
	//can be asked here rather than guessed.
	const CGLPixelFormatAttribute forcedSoftware[] = {
		kCGLPFAOpenGLProfile, static_cast< CGLPixelFormatAttribute >( kCGLOGLPVersion_GL4_Core ),
		kCGLPFARendererID, static_cast< CGLPixelFormatAttribute >( kCGLRendererGenericFloatID ),
		kCGLPFAColorSize, static_cast< CGLPixelFormatAttribute >( 24 ),
		kCGLPFAAlphaSize, static_cast< CGLPixelFormatAttribute >( 8 ),
		static_cast< CGLPixelFormatAttribute >( 0 )
	};
	CGLPixelFormatObj format = nullptr;
	GLint count              = 0;
	const char* renderer     = std::getenv( "CT_RENDERER" );
	if( renderer && std::string( renderer ) == "software" )
	{
		if( CGLChoosePixelFormat( forcedSoftware, &format, &count ) != kCGLNoError || format == nullptr )
			return nullptr;
	}
	else if( CGLChoosePixelFormat( accelerated, &format, &count ) != kCGLNoError || format == nullptr )
		if( CGLChoosePixelFormat( software, &format, &count ) != kCGLNoError || format == nullptr )
			return nullptr;
	CGLContextObj context = nullptr;
	const CGLError error  = CGLCreateContext( format, nullptr, &context );
	CGLDestroyPixelFormat( format );
	if( error != kCGLNoError )
		return nullptr;
	CGLSetCurrentContext( context );
	return context;
}

GLuint MakeTexture( int width, int height, const float* pixels )
{
	GLuint texture = 0;
	glGenTextures( 1, &texture );
	glBindTexture( GL_TEXTURE_2D, texture );
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, pixels );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glBindTexture( GL_TEXTURE_2D, 0 );
	return texture;
}

Floats ReadTexture( GLuint texture, int width, int height, GLenum format, int channels )
{
	Floats data( static_cast< size_t >( width ) * height * channels );
	glBindTexture( GL_TEXTURE_2D, texture );
	glPixelStorei( GL_PACK_ALIGNMENT, 1 );
	glGetTexImage( GL_TEXTURE_2D, 0, format, GL_FLOAT, data.data() );
	glBindTexture( GL_TEXTURE_2D, 0 );
	return data;
}

//---------------------------------------------------------------------------
namespace
{
uint32_t lowbias32( uint32_t x )
{
	x ^= x >> 16;
	x *= 0x7feb352dU;
	x ^= x >> 15;
	x *= 0x846ca68bU;
	x ^= x >> 16;
	return x;
}

double hash01( uint32_t a, uint32_t b = 0 )
{
	return static_cast< double >( lowbias32( a ^ lowbias32( b + 0x9e3779b9U ) ) ) / 4294967296.0;
}

void putU32( std::vector< unsigned char >& out, uint32_t value )
{
	out.push_back( static_cast< unsigned char >( value >> 24 ) );
	out.push_back( static_cast< unsigned char >( value >> 16 ) );
	out.push_back( static_cast< unsigned char >( value >> 8 ) );
	out.push_back( static_cast< unsigned char >( value ) );
}

void putChunk( std::vector< unsigned char >& out, const char* type, const std::vector< unsigned char >& data )
{
	putU32( out, static_cast< uint32_t >( data.size() ) );
	const size_t start = out.size();
	out.insert( out.end(), type, type + 4 );
	out.insert( out.end(), data.begin(), data.end() );
	uLong crc = crc32( 0L, Z_NULL, 0 );
	crc       = crc32( crc, out.data() + start, static_cast< uInt >( 4 + data.size() ) );
	putU32( out, static_cast< uint32_t >( crc ) );
}
} // namespace

/// A stained-glass card: a few dozen coloured cells with dark leading, a
/// sky-to-ground gradient under them, and fine hatching -- colour and detail
/// at every scale, so a plasma made of it shows where it has been stirred.
Floats BuildCard( int width, int height )
{
	Floats card( static_cast< size_t >( width ) * height * 4 );
	struct Site
	{
		double u, v, r, g, b;
	};
	std::vector< Site > sites;
	const double aspect = static_cast< double >( width ) / height;
	for( uint32_t i = 0; i < 48; ++i )
	{
		Site s;
		s.u = hash01( i, 1 ) * aspect;
		s.v = hash01( i, 2 );
		//Saturated, bright hues round the wheel.
		const double hue = hash01( i, 3 ) * 6.0;
		const double f   = hue - std::floor( hue );
		const int sector = static_cast< int >( hue ) % 6;
		const double on = 1.0, off = 0.15, up = off + ( on - off ) * f, down = on - ( on - off ) * f;
		const double rgb[ 6 ][ 3 ] = { { on, up, off }, { down, on, off }, { off, on, up },
			                           { off, down, on }, { up, off, on }, { on, off, down } };
		s.r = rgb[ sector ][ 0 ];
		s.g = rgb[ sector ][ 1 ];
		s.b = rgb[ sector ][ 2 ];
		sites.push_back( s );
	}
	for( int y = 0; y < height; ++y )
		for( int x = 0; x < width; ++x )
		{
			const double u = ( x + 0.5 ) / height, v = ( y + 0.5 ) / height;
			double best = 1e9, second = 1e9;
			const Site* nearest = &sites[ 0 ];
			for( const Site& s : sites )
			{
				const double d = ( u - s.u ) * ( u - s.u ) + ( v - s.v ) * ( v - s.v );
				if( d < best )
				{
					second  = best;
					best    = d;
					nearest = &s;
				}
				else if( d < second )
					second = d;
			}
			const double edge  = std::sqrt( second ) - std::sqrt( best );
			const double lead  = std::clamp( edge / 0.006, 0.0, 1.0 );
			const double hatch = 0.9 + 0.1 * std::sin( 400.0 * ( u + v ) );
			const size_t o     = ( static_cast< size_t >( y ) * width + x ) * 4;
			card[ o + 0 ]      = static_cast< float >( nearest->r * lead * hatch );
			card[ o + 1 ]      = static_cast< float >( nearest->g * lead * hatch );
			card[ o + 2 ]      = static_cast< float >( nearest->b * lead * hatch );
			card[ o + 3 ]      = 1.0f;
		}
	return card;
}

bool WritePng( const std::string& path, int width, int height, const Floats& rgba )
{
	std::vector< unsigned char > raw;
	raw.reserve( static_cast< size_t >( height ) * ( 1 + static_cast< size_t >( width ) * 4 ) );
	for( int y = height - 1; y >= 0; --y )
	{
		raw.push_back( 0 );
		for( int x = 0; x < width; ++x )
			for( int c = 0; c < 4; ++c )
			{
				const float v = rgba[ ( static_cast< size_t >( y ) * width + x ) * 4 + c ];
				raw.push_back( static_cast< unsigned char >( std::lround( std::clamp( v, 0.0f, 1.0f ) * 255.0f ) ) );
			}
	}
	uLongf size = compressBound( static_cast< uLong >( raw.size() ) );
	std::vector< unsigned char > compressed( size );
	if( compress2( compressed.data(), &size, raw.data(), static_cast< uLong >( raw.size() ), 6 ) != Z_OK )
		return false;
	compressed.resize( size );

	std::vector< unsigned char > png = { 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' };
	std::vector< unsigned char > ihdr;
	putU32( ihdr, static_cast< uint32_t >( width ) );
	putU32( ihdr, static_cast< uint32_t >( height ) );
	ihdr.insert( ihdr.end(), { 8, 6, 0, 0, 0 } );
	putChunk( png, "IHDR", ihdr );
	putChunk( png, "IDAT", compressed );
	putChunk( png, "IEND", {} );
	FILE* file = fopen( path.c_str(), "wb" );
	if( file == nullptr )
		return false;
	const size_t written = fwrite( png.data(), 1, png.size(), file );
	fclose( file );
	return written == png.size();
}

//---------------------------------------------------------------------------
Rig::~Rig()
{
	plugin.DeInitGL();
	if( outputFBO )
		glDeleteFramebuffers( 1, &outputFBO );
	if( outputTexture )
		glDeleteTextures( 1, &outputTexture );
	if( sourceTexture )
		glDeleteTextures( 1, &sourceTexture );
}

bool Rig::Init( int w, int h, const Floats* picture )
{
	//Under --size the output raster changes and the grid does not: the grid
	//keeps the aspect of the raster the check asked for, so every physics
	//check runs on exactly the grid it was written for, and only the light's
	//path to the raster -- emission, glare, composite -- sees the new size.
	if( g_rasterW > 0 && g_rasterH > 0 )
		plugin.SetGridRasterForTest( w, h );
	ChooseRaster( w, h );
	if( picture && picture->size() != static_cast< size_t >( w ) * h * 4 )
	{
		std::fprintf( stderr, "a check's picture is not %dx%d -- it must size it with ChooseRaster()\n", w, h );
		return false;
	}
	width  = w;
	height = h;
	FFGLViewportStruct viewport = {};
	viewport.width              = static_cast< FFUInt32 >( width );
	viewport.height             = static_cast< FFUInt32 >( height );
	if( plugin.InitGL( &viewport ) != FF_SUCCESS )
	{
		std::fprintf( stderr, "InitGL failed -- see the diagnostics log for which shader\n" );
		return false;
	}
	plugin.SetClockScaleForTest( 1.0 );

	const Floats card = picture ? *picture : BuildCard( width, height );
	sourceTexture     = MakeTexture( width, height, card.data() );
	outputTexture     = MakeTexture( width, height, nullptr );
	glGenFramebuffers( 1, &outputFBO );
	glBindFramebuffer( GL_FRAMEBUFFER, outputFBO );
	glFramebufferTexture2D( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, outputTexture, 0 );
	if( glCheckFramebufferStatus( GL_FRAMEBUFFER ) != GL_FRAMEBUFFER_COMPLETE )
	{
		std::fprintf( stderr, "the harness's own output framebuffer is not complete\n" );
		return false;
	}
	inputStruct.Width = inputStruct.HardwareWidth = static_cast< FFUInt32 >( width );
	inputStruct.Height = inputStruct.HardwareHeight = static_cast< FFUInt32 >( height );
	inputStruct.Handle                              = sourceTexture;
	inputs[ 0 ]                                     = &inputStruct;
	process.numInputTextures = 1;
	process.inputTextures    = inputs;
	process.HostFBO          = outputFBO;
	return true;
}

void Rig::Upload( const Floats& picture )
{
	glBindTexture( GL_TEXTURE_2D, sourceTexture );
	glTexSubImage2D( GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT, picture.data() );
	glBindTexture( GL_TEXTURE_2D, 0 );
}

void Rig::Set( unsigned int id, float value )
{
	plugin.SetFloatParameter( id, value );
}

void Rig::Press( unsigned int id )
{
	plugin.SetFloatParameter( id, 1.0f );
	plugin.SetFloatParameter( id, 0.0f );
}

void Rig::Quiet()
{
	Set( PT_FEED, 0.0f );
	Set( PT_CLIP_HEATS, 0.0f );
	Set( PT_FUEL, 0.0f );
	Set( PT_DRIVE, 0.0f );
	Set( PT_COIL_SPIN, 0.5f );
	Set( PT_CURVATURE, 0.0f );
	Set( PT_AUDIO_HEAT, 0.0f );
	Set( PT_AUDIO_PELLETS, 0.0f );
	Set( PT_RESISTIVITY, 0.0f );
	Set( PT_COOLING, 0.0f );
}

bool Rig::Render( int frames )
{
	for( int i = 0; i < frames; ++i )
	{
		const double seconds = static_cast< double >( frame ) / fps;
		plugin.SetTime( seconds );
		const double beat  = std::fmod( seconds, 0.5 );
		const float strike = feed == AudioFeed::Pulses ? static_cast< float >( 0.15 + 1.5 * std::exp( -beat / 0.06 ) ) : 0.0f;
		for( int bin = 0; bin < audio::kBins; ++bin )
		{
			const float across = static_cast< float >( bin ) / static_cast< float >( audio::kBins - 1 );
			const float shape  = 0.7f * ( 1.0f - across ) * ( 1.0f - across ) + 0.2f * ( 0.5f + 0.5f * std::sin( 25.0f * across ) );
			plugin.SetParamElementValue( PT_AUDIO, static_cast< unsigned int >( bin ), shape * strike );
		}
		++frame;

		glBindFramebuffer( GL_FRAMEBUFFER, outputFBO );
		glViewport( 0, 0, width, height );
		glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
		glClear( GL_COLOR_BUFFER_BIT );
		if( plugin.ProcessOpenGL( &process ) != FF_SUCCESS )
		{
			std::fprintf( stderr, "ProcessOpenGL failed\n" );
			return false;
		}
	}
	return true;
}

Floats Rig::Output() const
{
	Floats pixels( static_cast< size_t >( width ) * height * 4 );
	glBindFramebuffer( GL_FRAMEBUFFER, outputFBO );
	glPixelStorei( GL_PACK_ALIGNMENT, 1 );
	glReadPixels( 0, 0, width, height, GL_RGBA, GL_FLOAT, pixels.data() );
	return pixels;
}

//---------------------------------------------------------------------------
float ParamGeometric( double value, double low, double high )
{
	return static_cast< float >( std::log( value / low ) / std::log( high / low ) );
}
float DetailParam( int cells )
{
	for( int i = 0; i < kDetailCount; ++i )
		if( kDetailCells[ i ] == cells )
			return static_cast< float >( i );
	return 1.0f;
}
float PolesParam( int poles )
{
	for( int i = 0; i < kPoleOptions; ++i )
		if( kPoleCounts[ i ] == poles )
			return static_cast< float >( i );
	return 0.0f;
}
float CurvatureParam( double g )
{
	return static_cast< float >( std::sqrt( std::max( g, 0.0 ) / 3.0 ) );
}
float FieldParam( double field )
{
	return ParamGeometric( field, 0.25, 4.0 );
}
float GuideParam( double multiple )
{
	return static_cast< float >( multiple / 2.0 );
}
float TemperatureParam( double beta )
{
	return ParamGeometric( beta, 0.1, 50.0 );
}
float ResistivityParam( double eta )
{
	return eta <= 0.0 ? 0.0f : ParamGeometric( eta, 1e-5, 1e-2 );
}
float BallSizeParam( double radius )
{
	return static_cast< float >( ( radius - 0.04 ) / 0.36 );
}

//---------------------------------------------------------------------------
Snapshot Snapshot::Take( const ContainmentPlugin& plugin, double gamma )
{
	Snapshot s;
	s.nx    = plugin.CurrentGrid().nx;
	s.ny    = plugin.CurrentGrid().ny;
	s.dx    = plugin.CurrentGrid().dx;
	s.gamma = gamma;
	s.a     = ReadTexture( plugin.StateTextureID( 0 ), s.nx, s.ny );
	s.b     = ReadTexture( plugin.StateTextureID( 1 ), s.nx, s.ny );
	s.c     = ReadTexture( plugin.StateTextureID( 2 ), s.nx, s.ny );
	s.d     = ReadTexture( plugin.StateTextureID( 3 ), s.nx, s.ny );
	return s;
}

double Snapshot::P( int i, int j ) const
{
	const size_t o   = At( i, j );
	const double rho = a[ o ];
	const double kin = 0.5 * ( a[ o + 1 ] * a[ o + 1 ] + a[ o + 2 ] * a[ o + 2 ] + a[ o + 3 ] * a[ o + 3 ] ) / rho;
	const double mag = 0.5 * ( b[ o + 1 ] * b[ o + 1 ] + b[ o + 2 ] * b[ o + 2 ] + b[ o + 3 ] * b[ o + 3 ] );
	return ( gamma - 1.0 ) * ( b[ o ] - kin - mag );
}

bool Snapshot::Finite() const
{
	for( const Floats* v : { &a, &b, &c, &d } )
		for( float x : *v )
			if( !std::isfinite( x ) )
				return false;
	return true;
}

double Snapshot::Sum( int texture, int channel ) const
{
	const Floats* v[ 4 ] = { &a, &b, &c, &d };
	double sum           = 0.0;
	for( size_t i = static_cast< size_t >( channel ); i < v[ texture ]->size(); i += 4 )
		sum += ( *v[ texture ] )[ i ];
	return sum;
}

StateBuilder::StateBuilder( int x, int y, double g ) : nx( x ), ny( y ), gamma( g )
{
	const size_t n = static_cast< size_t >( nx ) * ny * 4;
	a.assign( n, 0.0f );
	b.assign( n, 0.0f );
	c.assign( n, 0.0f );
	d.assign( n, 0.0f );
}

void StateBuilder::Set( int i, int j, double rho, double u, double v, double w, double p, double bx, double by,
                        double bz, double t0, double t1, double t2, double t3 )
{
	const size_t o = ( static_cast< size_t >( j ) * nx + i ) * 4;
	a[ o + 0 ]     = static_cast< float >( rho );
	a[ o + 1 ]     = static_cast< float >( rho * u );
	a[ o + 2 ]     = static_cast< float >( rho * v );
	a[ o + 3 ]     = static_cast< float >( rho * w );
	b[ o + 0 ] = static_cast< float >( p / ( gamma - 1.0 ) + 0.5 * rho * ( u * u + v * v + w * w )
	                                   + 0.5 * ( bx * bx + by * by + bz * bz ) );
	b[ o + 1 ] = static_cast< float >( bx );
	b[ o + 2 ] = static_cast< float >( by );
	b[ o + 3 ] = static_cast< float >( bz );
	d[ o + 0 ] = static_cast< float >( rho * t0 );
	d[ o + 1 ] = static_cast< float >( rho * t1 );
	d[ o + 2 ] = static_cast< float >( rho * t2 );
	d[ o + 3 ] = static_cast< float >( rho * t3 );
}

} // namespace cttest
