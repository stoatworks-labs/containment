#pragma once

#include "Audio.h"
#include "Controls.h"
#include "PassBuffer.h"
#include "Physics.h"
#include "Shaders.h"
#include "StateBuffer.h"

#include <FFGLSDK.h>

// After FFGLSDK.h, which is where FFUInt32 comes from.
#include "StoatworksAboutParams.h"

#include <string>
#include <vector>

namespace containment
{
/**
    Deliberate departures from the shipped model, for the harness's negative
    controls and for the tests that need a textbook setup (Brio-Wu's gamma 2,
    the Alfven wave's periodic box, a planar slab's far-away gravity centre).
    A default-constructed one is the shipped model exactly.
*/
struct TestModel
{
	float gamma         = kGamma;
	int boundary        = -1;  ///< -1: the Boundary parameter; else a Boundary value (2 = periodic)
	int solver          = 0;   ///< 0 HLLD, 1 HLL
	bool glm            = true;
	bool floors         = true;
	bool entropy        = true; ///< the dual-energy switch
	bool gravityCentre  = false;///< use the centre below instead of the ball's
	double gravityX     = 0.0;
	double gravityY     = 0.0;
	double gravityCore  = -1.0; ///< < 0: the shipped taper
	bool uniformGravity = false;///< g_eff not scaled by T / T_ref (the textbook RT slab)
	bool additiveGlow   = false;///< the glare ADDED rather than moved (the wrong model)
	bool legacyOpen     = false;///< Open as it was: no margin, no absorbing layer (the wrong model)
	float backgroundDensity  = kBackgroundDensity;
	float backgroundPressure = kBackgroundPressure;
};

/**
    The plugin.

    One model -- 2.5-D compressible ideal MHD on the GPU -- and one light (the
    plasma's own emission). See AGENTS.md for why each is built as it is.
*/
class ContainmentPlugin : public CFFGLPlugin
{
public:
	ContainmentPlugin();

	FFResult InitGL( const FFGLViewportStruct* viewport ) override;
	FFResult ProcessOpenGL( ProcessOpenGLStruct* input ) override;
	FFResult DeInitGL() override;

	FFResult SetFloatParameter( unsigned int index, float value ) override;
	float GetFloatParameter( unsigned int index ) override;
	char* GetTextParameter( unsigned int index ) override;
	/// The base class fails, and a failed default deletes the instance: without
	/// this no real host can load the plugin. See StoatworksAboutParams.h.
	FFResult SetTextParameter( unsigned int index, const char* value ) override;
	FFResult SetTime( double time ) override;

	//-------------------------------------------------------------------
	// For the harness. Nothing in the plugin's own operation calls these.
	//-------------------------------------------------------------------

	/// Replace one exact substring of every shader before the next InitGL: the
	/// one-character mutation that proves the harness drives shipped code.
	static void SetShaderMutationForTest( const std::string& find, const std::string& replace );

	void SetClockScaleForTest( double scale )
	{
		clockScale = scale;
	}
	void SetModelForTest( const TestModel& model )
	{
		test = model;
	}

	/// A parameter's effective value: the operator's, or the preset's where
	/// the Preset dropdown is on anything but Custom.
	float P( unsigned int index ) const;
	float EffectiveForTest( unsigned int index ) const
	{
		return P( index );
	}

	const Grid& CurrentGrid() const
	{
		return grid;
	}
	double SimTime() const
	{
		return simTime;
	}
	/// The coils as of the last frame.
	const Coils& CurrentCoils() const
	{
		return coils;
	}

	/// The state texture `index` (0..3, A..D) as it stands.
	GLuint StateTextureID( int index ) const
	{
		return state[ current ].Texture( index );
	}
	/// Overwrite the state, four RGBA float arrays of grid size, row 0 first.
	/// Signal speeds are recomputed on the next step.
	void LoadStateForTest( const std::vector< float >& a, const std::vector< float >& b, const std::vector< float >& c,
	                       const std::vector< float >& d );
	/// Advance the plasma by `duration` of simulated time with the current
	/// parameters and no events, planning substeps from a synchronous
	/// readback. Returns the substeps taken, or -1 if the clock stopped (a
	/// non-finite state, or `limit` substeps without getting there).
	int StepForTest( double duration, int limit = 200000 );

	int EmissionWidth() const
	{
		return grid.fx;
	}
	int EmissionHeight() const
	{
		return grid.fy;
	}
	GLuint EmissionTextureID() const
	{
		return emission.TextureID();
	}
	GLuint GlowTextureID( int stage ) const
	{
		return glow[ stage ].TextureID();
	}
	int GlowWidth() const
	{
		return glowWidth;
	}
	int GlowHeight() const
	{
		return glowHeight;
	}
	float GlowFraction() const
	{
		return glowFraction;
	}
	GLuint PotentialTextureID() const
	{
		return levels.empty() ? 0 : levels[ 0 ].solution[ levels[ 0 ].current ].TextureID();
	}

	/// Totals since the plugin started.
	long long SubstepsTaken() const
	{
		return substepsTaken;
	}
	/// Floor clamps, summed over cells and substeps, as read back.
	double FloorHits() const
	{
		return floorHits;
	}
	/// Cells that took their pressure from the entropy (the dual-energy
	/// switch), summed over cells and substeps, as read back.
	double EntropyCells() const
	{
		return entropyCells;
	}
	/// Frames in which the substep cap stopped simulated time keeping up.
	int CappedFrames() const
	{
		return cappedFrames;
	}

private:
	/// Everything the passes need, from the parameters.
	struct Model
	{
		double field = 1.0, guide = 0.0, spin = 0.0, curvature = 0.0;
		double eta = 0.0, cooling = 0.0, drive = 0.0, driveScale = 0.2;
		int poles = 0, boundary = 0;
		double coilRadius = 1.3;
		double ballX = 0.5, ballY = 0.5, ballRadius = 0.2, pressure = 0.25;
	};
	Model CurrentModel() const;

	void UpdateClock();
	bool EnsureBuffers( int width, int height, const Grid& wanted );

	/// Uniforms every state pass shares: grid, gas, boundary, coils, forces.
	void SetStateUniforms( GLuint program ) const;
	void DrawInto( const StateBuffer& target );

	void Ignite( GLuint input, float maxU, float maxV );
	void Sources( GLuint input, float maxU, float maxV, double hostDt, double simDt, int pellets, bool always = false );
	/// Reduce the state's speeds into the clock for substep k of n.
	void Clock( int substep, int substeps, double target );
	void Substep();
	/// Plan and run the steps for `target` of simulated time.
	void Advance( double target, bool synchronous );
	/// Read back the last clock: ( dt, done, c_h, floors ).
	void ReadClock( float out[ 4 ] );

	void Emission();
	void Glow();
	void Potential();

	float params[ PT_COUNT ] = {};
	TestModel test;

	ffglex::FFGLShader programs[ static_cast< int >( Program::Count ) ];
	ffglex::FFGLScreenQuad quad;
	GLuint Id( Program p ) const
	{
		return programs[ static_cast< int >( p ) ].GetGLID();
	}

	Grid grid;
	StateBuffer state[ 2 ];
	StateBuffer star;
	StateBuffer fluxX;///< ( nx + 1 ) x ny: the x faces
	StateBuffer fluxY;///< nx x ( ny + 1 ): the y faces
	int current = 0;
	std::vector< PassBuffer > reduction;
	PassBuffer clock[ 2 ];
	int clockIndex = 0;

	PassBuffer emission;
	std::vector< PassBuffer > glowReduce;
	PassBuffer glowTemp;
	PassBuffer glow[ kGlowStages ];
	int glowWidth = 0, glowHeight = 0;
	float glowFraction = 0.0f;

	struct Level
	{
		PassBuffer solution[ 2 ];
		PassBuffer rhs;
		PassBuffer residual;
		int current = 0;
		int nx = 0, ny = 0;
		double spacing = 0.0;
	};
	std::vector< Level > levels;
	bool havePotential = false;

	Coils coils;
	Coils previousCoils;///< as of the last frame: the wall is driven by the change
	bool carryCoils = false;///< this frame's change of the coils is carried through the vessel (Wall)
	Drive drive;
	double spinAngle     = 0.0;
	double coilStrength  = 1.0;
	int frameSubsteps    = 0;
	float lastSpeed      = 0.0f;///< max summed signal speed, as last read back
	bool speedKnown      = false;
	bool clockPending    = false;///< a frame's clock waits to be read back
	double lastTarget    = 0.0;

	long long substepsTaken = 0;
	double floorHits        = 0.0;
	double entropyCells     = 0.0;
	float lastEntropyCells  = 0.0f;
	int cappedFrames        = 0;
	int framesSinceCapLog   = 1 << 20;

	//-------------------------------------------------------------------
	// Time. See rosette: the host's clock unit is voted on against the wall
	// clock, because Resolume has sent both seconds and milliseconds.
	//-------------------------------------------------------------------
	double hostTime     = -1.0;
	double lastRawTime  = -1.0;
	double lastWallTime = -1.0;
	double wallStart    = -1.0;
	double clockScale   = 0.0;
	int secondsVotes    = 0;
	int millisVotes     = 0;
	int clockFrames     = 0;
	double now          = 0.0;
	double lastNow      = -1.0;
	double simTime      = 0.0;
	bool settledJump    = false;

	//-------------------------------------------------------------------
	// Events, on the rising edge.
	//-------------------------------------------------------------------
	bool ignitePending = true;///< the first frame lays down a ball
	int pelletPresses  = 0;
	bool igniteHeld    = false;
	bool pelletHeld    = false;

	audio::Analyser analyser;
};

} // namespace containment
