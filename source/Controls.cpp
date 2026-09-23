#include "Controls.h"

#include <algorithm>
#include <cmath>

namespace containment
{
namespace
{
float clamp01( float value )
{
	return std::clamp( value, 0.0f, 1.0f );
}

/// value^0 = low, value^1 = high, with every octave the same width of travel.
float geometric( float value, float low, float high )
{
	return low * std::pow( high / low, clamp01( value ) );
}

float linear( float value, float low, float high )
{
	return low + ( high - low ) * clamp01( value );
}

constexpr float kSpeedLow  = 0.02f;
constexpr float kSpeedHigh = 2.0f;
} // namespace

float BallSizeFromParam( float value )
{
	return linear( value, 0.04f, 0.4f );
}

float TemperatureFromParam( float value )
{
	return geometric( value, 0.1f, 50.0f );
}

float FeedFromParam( float value )
{
	return clamp01( value );
}

float ClipHeatsFromParam( float value )
{
	return clamp01( value );
}

float FuelFromParam( float value )
{
	return linear( value, 0.0f, 2.0f );
}

float FieldFromParam( float value )
{
	return geometric( value, 0.25f, 4.0f );
}

float GuideFieldFromParam( float value )
{
	return linear( value, 0.0f, 2.0f );
}

float CoilRadiusFromParam( float value )
{
	return linear( value, 1.1f, 3.0f );
}

float CoilSpinFromParam( float value )
{
	return linear( value, -1.0f, 1.0f );
}

float CurvatureFromParam( float value )
{
	const float v = clamp01( value );
	return 3.0f * v * v;
}

float SpeedFromParam( float value )
{
	//Exactly zero at the bottom: a frozen plasma is a thing an operator wants,
	//and "very slow" is not the same thing.
	if( value <= 0.0f )
		return 0.0f;
	return geometric( value, kSpeedLow, kSpeedHigh );
}

float ParamFromSpeed( float speed )
{
	if( speed <= 0.0f )
		return 0.0f;
	const float clamped = std::clamp( speed, kSpeedLow, kSpeedHigh );
	return std::log( clamped / kSpeedLow ) / std::log( kSpeedHigh / kSpeedLow );
}

float ResistivityFromParam( float value )
{
	if( value <= 0.0f )
		return 0.0f;
	return geometric( value, 1.0e-5f, 1.0e-2f );
}

float CoolingFromParam( float value )
{
	const float v = clamp01( value );
	return 8.0f * v * v * v;
}

float DriveFromParam( float value )
{
	const float v = clamp01( value );
	return 4.0f * v * v;
}

float DriveScaleFromParam( float value )
{
	return geometric( value, 0.05f, 0.6f );
}

float AudioHeatFromParam( float value )
{
	return linear( value, 0.0f, 4.0f );
}

float ExposureFromParam( float value )
{
	return linear( value, -4.0f, 8.0f );
}

float GlowFromParam( float value )
{
	return linear( value, 0.0f, 0.9f );
}

int LineCountFromParam( float value )
{
	return std::clamp( static_cast< int >( std::lround( linear( value, 4.0f, 48.0f ) ) ), 4, 48 );
}

} // namespace containment
