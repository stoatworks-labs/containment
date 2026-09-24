#pragma once

#include "Harness.h"

#include <functional>

namespace cttest
{
// Each returns 0 on pass; the failures are counted in g_failures.
int RunBrioWu( const Perturb& perturb );
int RunAlfven( const Perturb& perturb );
int RunConserve( const Perturb& perturb );
int RunDivB( const Perturb& perturb );
int RunBalance( const Perturb& perturb );
int RunRT( const Perturb& perturb );
int RunCusp( const Perturb& perturb );
int RunFrozen( const Perturb& perturb );
int RunQuench( const Perturb& perturb );
int RunResist( const Perturb& perturb );
int RunFloors( const Perturb& perturb );
int RunStill( const Perturb& perturb );
int RunGlow( const Perturb& perturb );
int RunState( const Perturb& perturb );
int RunMutation( const Perturb& perturb );
int RunEquilibrium( const Perturb& perturb );
int RunOpen( const Perturb& perturb );
int RunPresets( const Perturb& perturb );
int RunReference( const Perturb& perturb );
int RunVacuum( const Perturb& perturb );
int RunNames( const Perturb& perturb );
/// The negative controls; `offline` runs only those of the checks above that
/// need no GL context.
int RunNegative( bool offline = false );
/// `before` is applied to every rig first (main's --set / --preset).
int RunBench( const std::function< void( Rig& ) >& before = {} );
} // namespace cttest
