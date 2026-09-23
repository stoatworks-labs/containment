#pragma once

#include "Harness.h"

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
int RunNegative();
int RunBench();
} // namespace cttest
