#pragma once

/**
    The host's parameters, and what they mean in physical units.

    Every numeric parameter the host sees is a plain 0..1 float, because
    `SetParamInfo` clamps an `FF_TYPE_STANDARD` default into 0..1 *before*
    `SetParamRange` could widen it. The conversions all live in Controls.cpp,
    one function per control, and the shaders are handed the physical value.
    Option, boolean and event parameters hold the element value itself.

    ------------------------------------------------------------------ units

    The plasma runs in the units ideal MHD is usually written in:

    - **length**: the frame height is 1. The grid is square-celled, Detail
      cells on the short side, so every physics check is raster-independent
      by construction;
    - **mu0 = 1**, so magnetic pressure is B^2/2 and the Alfven speed B/sqrt(rho);
    - **the reference state** is rho0 = 1 and B_ref = 1, so the unit of time
      is one Alfven crossing of the frame height at the reference field and
      density, tau_A = L sqrt(rho0) / B_ref = 1;
    - **gamma = 5/3**, a monatomic plasma.

    Field is the bottle's strength in units of B_ref, and Temperature is the
    ball's beta against B_ref (p0 = beta0 B_ref^2 / 2). They are separate on
    purpose: turning Field up confines the SAME ball harder, rather than
    rescaling the whole picture.
*/

namespace containment
{
/**
    Parameter ids.

    **Append only.** `SetParamGroup` collapses runs of consecutive same-group
    ids, so inserting an id mid-enum silently splits a group in two; and every
    saved composition stores parameters by index.
*/
/// What a host shows. The FFGL name field is char[ 16 ] and NOT null-
/// terminated, so a longer name is truncated without a word; `cttest --names`
/// holds this to 16 and `oxbow probe` reads the bundle's copy back.
constexpr const char* kDisplayName = "SW Containment";
constexpr const char* kPluginCode  = "CT01";

enum ParamId : unsigned int
{
	// -- Preset -------------------------------------------------------------
	// Element 0 is Custom; any other lays a row of Presets.h over the
	// operator's values at read time (the fleet's override model).
	PT_PRESET = 0,

	// -- Ball ---------------------------------------------------------------
	PT_IGNITE,
	PT_BALL_SIZE,
	PT_BALL_X,
	PT_BALL_Y,
	PT_TEMPERATURE,
	PT_PROFILE,
	PT_PELLET,
	PT_FEED,
	PT_CLIP_HEATS,
	PT_FUEL,

	// -- Bottle -------------------------------------------------------------
	PT_FIELD,
	PT_GUIDE_FIELD,
	PT_POLES,
	PT_COIL_RADIUS,
	PT_COIL_SPIN,
	PT_CURVATURE,
	PT_QUENCH,
	PT_BOUNDARY,

	// -- Plasma -------------------------------------------------------------
	PT_SPEED,
	PT_RESISTIVITY,
	PT_COOLING,
	PT_DETAIL,
	PT_DRIVE,
	PT_DRIVE_SCALE,

	// -- Audio --------------------------------------------------------------
	PT_AUDIO,
	PT_AUDIO_HEAT,
	PT_AUDIO_PELLETS,

	// -- Light --------------------------------------------------------------
	PT_EXPOSURE,
	PT_TINT,
	PT_RAMP,
	PT_GLOW,
	PT_FIELD_LINES,
	PT_LINE_COUNT,
	PT_VIEW,
	PT_MIX,

	// -- The Stoatworks About block -----------------------------------------
	// One display-only text line, then one button per link the block carries.
	// Containment.cpp static_asserts this run against `about::kParamCount`.
	PT_ABOUT_TEXT,
	PT_ABOUT_BUTTON_1,
	PT_ABOUT_BUTTON_2,
	PT_ABOUT_BUTTON_3,
	PT_ABOUT_BUTTON_4,
	PT_COUNT
};

enum class Profile
{
	Gaussian = 0,
	TopHat,
	Count
};

enum class Boundary
{
	Open = 0,///< zero-gradient outflow: plasma that leaves the frame is gone; B is the coils'
	Wall,    ///< perfectly conducting and reflecting: nothing crosses it
	Count,
	Periodic = 2///< harness only: the Alfven wave needs a periodic box. Not a host option.
};

enum class Ramp
{
	Hot = 0,///< deep red -> violet -> white
	Aurora, ///< the 557.7 nm oxygen line, going white-green when hot
	Count
};

enum class View
{
	Picture = 0,
	Density,
	Pressure,
	FieldStrength,
	Beta,
	Speed,
	DivB,
	Count
};

/// The multipole orders Poles offers. 0 is the guide field alone.
constexpr int kPoleCounts[] = { 0, 4, 6, 8, 12 };
constexpr int kPoleOptions  = 5;

/// Grid cells on the frame's short side.
constexpr int kDetailCells[] = { 128, 256, 512, 1024 };
constexpr int kDetailCount   = 4;

//---------------------------------------------------------------------------
// Constants of the model. None of these is a control.
//---------------------------------------------------------------------------

/// The ratio of specific heats: a fully ionised, monatomic plasma.
constexpr float kGamma = 5.0f / 3.0f;

/// The plasma outside the ball. Tenuous and cold, so it emits almost nothing
/// (rho^2 sqrt(T) is 0.2% of the ball's), but not vacuum: an Alfven speed of
/// B / sqrt(rho) in real vacuum is infinite, and the CFL step with it.
constexpr float kBackgroundDensity  = 0.1f;
constexpr float kBackgroundPressure = 0.005f;

/// The floors. Vacuum cannot be represented, and an MHD pressure is a small
/// difference of large numbers wherever beta is small -- so p can come out
/// negative by rounding alone. These are the values it is clamped to, and
/// every clamp is counted (Diag, and `cttest --floors`). They are far below
/// anything the model is meant to reach: a floor that fires often means the
/// model is wrong, and is not a knob for the look.
constexpr float kDensityFloor  = 1.0e-4f;
constexpr float kPressureFloor = 1.0e-6f;

/// The Courant number, on the SUM of the two directions' signal speeds
/// (unsplit MUSCL-Hancock is stable for sum <= 1).
constexpr float kCourant = 0.8f;

/// The most substeps one frame may take. Past this, simulated time runs slow
/// -- never unstable. Logged through Diag when it bites.
constexpr int kMaxSubsteps = 48;

/// Open's margin: simulated cells all round the frame, as a fraction of the
/// frame height, holding the absorbing layer. Never shown.
constexpr float kMarginFraction = 0.1f;

/// The absorbing layer: e-foldings the fastest wave loses crossing it once.
constexpr float kSpongeEFolds = 6.0f;

/// The absorbing layer's rate rises as depth^kSpongePower: zero, with zero
/// slope, where the frame ends.
constexpr float kSpongePower = 2.0f;

/// The dual-energy switch: where the pressure the total energy gives is less
/// than this fraction of the kinetic plus magnetic energy, the pressure comes
/// from the entropy carried alongside instead. See the update pass.
constexpr float kEntropySwitch = 0.02f;

/// GLM damping, Mignone & Tzeferacos (2010): psi *= exp( -alpha c_h dt / dx ).
constexpr float kGLMAlpha = 0.4f;

/// The pellet: a cold, dense blob at the ball's centre.
constexpr float kPelletDensity = 2.0f;
constexpr float kPelletRadius  = 0.035f;

/// How long the coils take to lose their current after a Quench, in units of
/// tau_A. The boundary field decays as exp( -t / kQuenchTime ).
constexpr float kQuenchTime = 0.15f;

/// The stirring force's correlation time, in tau_A, and its mode count.
constexpr float kDriveTime  = 0.6f;
constexpr int kDriveModes   = 12;
/// The stirring's reach: a Gaussian window of this many ball radii about the
/// ball's centre (the stream function is windowed, so the force stays
/// divergence-free).
constexpr float kDriveWindow = 1.5f;

/// Emission rho^2 sqrt(T) at the reference ball (rho 1, T 1/4), which
/// Exposure 0 EV maps to 1.
constexpr float kEmissionReference = 0.5f;

//---------------------------------------------------------------------------
// The mappings. Each says its range and its shape.
//---------------------------------------------------------------------------

/// Ball radius, frame heights: 0.04 to 0.4, linearly. For the Gaussian
/// profile it is the 1/e radius; for the top hat, the edge.
float BallSizeFromParam( float value );

/// The ball's beta against B_ref, 0.1 to 50, geometrically; p0 = beta0 / 2.
float TemperatureFromParam( float value );

/// Feed as the fraction of the way the tracer moves towards the clip in one
/// frame at 60 fps. 0 to 1, linearly.
float FeedFromParam( float value );

/// Clip Heats: 0 to 1, linearly. At 1 the brightest pixel starts at twice the
/// ball's pressure and the darkest at none.
float ClipHeatsFromParam( float value );

/// Fuel: a steady gas puff that tops the ball's footprint back up towards its
/// ignition profile, at this rate per tau_A: 0 to 2, linearly. It only adds
/// (mass at rest, heat), never removes, so it holds a leaking ball up rather
/// than pinning it. At 0 the bottle runs down.
float FuelFromParam( float value );

/// The bottle's field in units of B_ref: 0.25 to 4, geometrically. For a
/// cusp it is |B| at radius 0.5 (the rim of the frame's inscribed circle).
float FieldFromParam( float value );

/// Guide field Bz as a multiple of Field: 0 to 2, linearly.
float GuideFieldFromParam( float value );

/// Coil radius, frame heights from the centre: 1.1 to 3, linearly. Never
/// inside the frame: the plugin pushes it out to 1.05x the half-diagonal.
float CoilRadiusFromParam( float value );

/// Coil rotation rate, radians per tau_A: -1 to 1, linearly, 0 at the middle.
float CoilSpinFromParam( float value );

/// Effective gravity g_eff at the ball's own temperature, B_ref^2 / ( rho0 L ):
/// 0 to 3, quadratically. Elsewhere it scales as T / T_ball (the curvature
/// drift goes as T / R_c), so its force density is Curvature p / T_ball.
float CurvatureFromParam( float value );

/// Speed: Alfven crossing times per second of host time. 0.02 to 2,
/// geometrically, and exactly 0 at the bottom of the travel.
float SpeedFromParam( float value );
float ParamFromSpeed( float speed );

/// Resistivity eta (magnetic diffusivity, mu0 = 1): 0, then 1e-5 to 1e-2,
/// geometrically.
float ResistivityFromParam( float value );

/// Bremsstrahlung cooling coefficient Lambda, dE/dt = -Lambda rho^2 sqrt(T):
/// 0 to 8, cubically.
float CoolingFromParam( float value );

/// Drive: the stirring acceleration's rms, B_ref^2 / ( rho0 L ): 0 to 4, quadratically.
float DriveFromParam( float value );

/// Drive Scale: the stirring wavelength, frame heights: 0.05 to 0.6, geometrically.
float DriveScaleFromParam( float value );

/// Audio Heat: heating power at full level, in p0 per tau_A: 0 to 4, linearly.
float AudioHeatFromParam( float value );

/// Exposure in EV: -4 to +8, linearly. 0 EV maps the reference ball to 1.
float ExposureFromParam( float value );

/// Glow: the fraction of the light the camera's glare moves out of the core,
/// 0 to 0.9, linearly.
float GlowFromParam( float value );

/// Field Line Count: 4 to 48 contours across the bottle's flux, linearly.
int LineCountFromParam( float value );

} // namespace containment
