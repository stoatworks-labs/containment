#include "StateBuffer.h"

namespace containment
{
bool StateBuffer::Ensure( int w, int h )
{
	if( w <= 0 || h <= 0 )
		return false;
	if( framebuffer != 0 && w == width && h == height )
		return true;

	Destroy();

	GLint previousFBO = 0, previousTexture = 0, previousViewport[ 4 ] = {};
	glGetIntegerv( GL_FRAMEBUFFER_BINDING, &previousFBO );
	glGetIntegerv( GL_TEXTURE_BINDING_2D, &previousTexture );
	glGetIntegerv( GL_VIEWPORT, previousViewport );

	glGenTextures( kTargets, textures );
	for( GLuint texture : textures )
	{
		glBindTexture( GL_TEXTURE_2D, texture );
		glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	}
	glBindTexture( GL_TEXTURE_2D, static_cast< GLuint >( previousTexture ) );

	glGenFramebuffers( 1, &framebuffer );
	glBindFramebuffer( GL_FRAMEBUFFER, framebuffer );
	const GLenum attachments[ kTargets ] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2,
		                                     GL_COLOR_ATTACHMENT3 };
	for( int i = 0; i < kTargets; ++i )
		glFramebufferTexture2D( GL_FRAMEBUFFER, attachments[ i ], GL_TEXTURE_2D, textures[ i ], 0 );
	//The draw-buffer list is framebuffer state, so this holds for every pass
	//that binds it.
	glDrawBuffers( kTargets, attachments );
	const bool complete = glCheckFramebufferStatus( GL_FRAMEBUFFER ) == GL_FRAMEBUFFER_COMPLETE;
	if( complete )
	{
		glViewport( 0, 0, w, h );
		glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
		glClear( GL_COLOR_BUFFER_BIT );
	}
	glBindFramebuffer( GL_FRAMEBUFFER, static_cast< GLuint >( previousFBO ) );
	glViewport( previousViewport[ 0 ], previousViewport[ 1 ], previousViewport[ 2 ], previousViewport[ 3 ] );

	if( !complete )
	{
		Destroy();
		return false;
	}
	width  = w;
	height = h;
	return true;
}

void StateBuffer::Destroy()
{
	if( framebuffer != 0 )
		glDeleteFramebuffers( 1, &framebuffer );
	framebuffer = 0;
	if( textures[ 0 ] != 0 )
		glDeleteTextures( kTargets, textures );
	for( GLuint& texture : textures )
		texture = 0;
	width = height = 0;
}

} // namespace containment
