#pragma once

/**
    Factory presets: a whole bottle an operator can reach in one gesture.

    The fleet's model (graticule's). **Presets are an OVERRIDE, not a write.**
    Resolume does not consume value events, so a plugin cannot push a preset's
    values back into the inspector; if it changed its own parameters the
    sliders would keep showing the old numbers. So while the Preset dropdown is
    on anything but Custom, the row's values are laid over the operator's at
    read time (`ContainmentPlugin::P`), and for those columns the inspector is
    not the truth. Element 0 of the dropdown is Custom and is not in this
    table: it means "the controls are the truth".

    **Standard parameters hold the host-facing 0..1; option and boolean
    parameters hold their real value** (Poles 2 is six poles, Boundary 0 is
    Open, Ramp 1 is Aurora). `cttest --presets` fails a row with a fraction in
    a discrete column, which would otherwise round without a word.

    Choosing a preset re-ignites: a new bottle wants a fresh ball.

    What a row leaves alone -- where the ball is, the Detail, Speed, the audio,
    View and Mix -- stays the operator's.

    ## The first row is also the constructor's defaults

    `Clip Orb` and the defaults in `Containment.cpp` are the same bottle,
    written twice, and `cttest --presets` fails when they drift apart.
*/

namespace containment
{
namespace presets
{
enum Param
{
	kBallSize,
	kTemperature,
	kProfile,
	kFeed,
	kClipHeats,
	kFuel,
	kField,
	kGuideField,
	kPoles,
	kCoilRadius,
	kCoilSpin,
	kCurvature,
	kQuench,
	kBoundary,
	kResistivity,
	kCooling,
	kDrive,
	kDriveScale,
	kExposure,
	kTint,
	kRamp,
	kGlow,
	kFieldLines,
	kLineCount,
	kParamCount
};

struct Preset
{
	const char* name;
	float v[ kParamCount ];
};

// clang-format off
//                                size   temp   prof feed  clip  fuel  field guide poles radius spin  curv   q  bnd res cool drive dscale expo    tint ramp glow   lines count
inline constexpr Preset kPresets[] = {
	{ "Clip Orb",             { 0.389f, 0.315f, 0, 0.1f, 0.0f, 0.15f, 0.25f, 0.7f, 2, 0.105f, 0.55f, 0.316f, 0, 0, 0, 0, 0.274f, 0.442f, 0.4167f, 0.0f, 0, 0.333f, 0.0f,  0.273f } },
	//The defaults: a ball held by the guide field and shaped by a weak
	//six-pole cusp, the coils turning, curvature and stirring just past
	//where the edge stays smooth, fuelled so the open bottle never runs
	//down. The clip keeps its own colours.

	{ "Green Orb",            { 0.389f, 0.315f, 0, 0.1f, 0.0f, 0.15f, 0.25f, 0.7f, 2, 0.105f, 0.62f, 0.36f,  0, 0, 0, 0, 0.33f,  0.4f,   0.45f,   0.9f, 1, 0.62f,  0.35f, 0.25f } },
	//The same bottle turning faster and stirred harder, in the oxygen green
	//of the Aurora ramp, with the glare and the field lines lit by the plasma
	//on them.

	{ "Guide Field Bubble",   { 0.306f, 0.4f,   1, 0.1f, 0.0f, 0.0f,  0.5f,  0.5f, 0, 0.105f, 0.5f,  0.0f,   0, 0, 0, 0, 0.0f,   0.442f, 0.4167f, 0.5f, 0, 0.4f,   0.0f,  0.273f } },
	//A theta-pinch cross-section: a top-hat ball in a uniform axial field,
	//no cusp, no curvature. It swells, rings and settles into the
	//diamagnetic bubble `cttest --balance` checks.

	{ "Cusp Leak",            { 0.222f, 0.482f, 0, 0.1f, 0.0f, 0.25f, 0.5f,  0.0f, 2, 0.105f, 0.5f,  0.0f,   0, 0, 0, 0, 0.0f,   0.442f, 0.4167f, 0.4f, 0, 0.45f,  0.3f,  0.273f } },
	//A hot ball in a six-pole cusp with no guide field: it squirts out along
	//the six field lines that lead out, at the angles `cttest --cusp` checks,
	//and the fuel keeps it squirting.

	{ "Rayleigh-Taylor",      { 0.389f, 0.315f, 1, 0.1f, 0.0f, 0.15f, 0.25f, 0.7f, 0, 0.105f, 0.5f,  0.7f,   0, 0, 0, 0, 0.1f,   0.442f, 0.4167f, 0.3f, 0, 0.4f,   0.0f,  0.273f } },
	//Strong curvature on a guide-field ball: nothing in the plane holds the
	//edge, so every wavelength is unstable and it breaks into fingers.

	{ "Quench Fireball",      { 0.306f, 0.741f, 1, 0.0f, 0.0f, 0.0f,  0.5f,  0.7f, 2, 0.105f, 0.5f,  0.0f,   1, 0, 0, 0, 0.0f,   0.442f, 0.4167f, 1.0f, 0, 0.6f,   0.0f,  0.273f } },
	//The coils off, a very hot ball: it free-expands into a fireball whose
	//front `cttest --quench` checks against the exact Riemann solution.
	//Press Ignite for another.
};
// clang-format on

inline constexpr int kCount = static_cast< int >( sizeof( kPresets ) / sizeof( kPresets[ 0 ] ) );

} // namespace presets
} // namespace containment
