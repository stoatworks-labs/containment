#pragma once

#include <FFGLSDK.h>

namespace containment
{
/**
    The plasma's state: four RGBA32F textures behind one framebuffer, so a pass
    can write all thirteen conserved quantities at once (multiple render
    targets). FFGLFBO has one colour attachment, so this is written against GL
    directly.

    32-bit everywhere, never half: this GPU truncates float -> half on store
    (millpond measured it), and a state fed back thousands of times a second
    would drift by half a step every step.

    Nearest filtering, clamp to edge: the solver reads texel for texel, and the
    boundary is applied in the shader, never by a wrap mode.
*/
class StateBuffer
{
public:
	static constexpr int kTargets = 4;

	/// Allocate at this size, reusing the buffer if it already matches, and
	/// clear new storage to zero. Binds nothing it does not put back.
	bool Ensure( int width, int height );
	void Destroy();

	GLuint Texture( int index ) const
	{
		return textures[ index ];
	}
	GLuint Framebuffer() const
	{
		return framebuffer;
	}
	int Width() const
	{
		return width;
	}
	int Height() const
	{
		return height;
	}

private:
	GLuint textures[ kTargets ] = {};
	GLuint framebuffer          = 0;
	int width                   = 0;
	int height                  = 0;
};

} // namespace containment
