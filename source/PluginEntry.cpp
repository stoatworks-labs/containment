#include "Containment.h"

/**
    The one registration.

    This file is listed directly in the Containment MODULE target, not in
    containment_core: `CFFGLPluginInfo` registers itself from a file-scope
    constructor and nothing ever references it by name, so in a STATIC archive
    the linker is entitled to drop the whole translation unit -- giving a
    bundle that loads, exports `plugMain`, and reports that it contains no
    plugins. The core stays an OBJECT library for the same reason.

    The name is `SW Containment`, fourteen characters. The FFGL name field is
    `char[ 16 ]` and is **not** null-terminated, so the host truncates without
    saying so. `oxbow probe` is what reads it back the way a host does.
*/
namespace
{
class ContainmentEffect : public containment::ContainmentPlugin
{
};
} // namespace

static CFFGLPluginInfo PluginInfo(
	PluginFactory< ContainmentEffect >,                   // Create method
	"CT01",                                               // Plugin unique ID of maximum length 4
	"SW Containment",                                     // Plugin name
	2,                                                    // API major version number
	1,                                                    // API minor version number
	0,                                                    // Plugin major version number
	1,                                                    // Plugin minor version number
	FF_EFFECT,                                            // Plugin type
	"A ball of hot plasma held in a magnetic bottle, obeying ideal magnetohydrodynamics. The clip is what the "
	"plasma is made of: it swells against the field, rings, leaks through the cusps and breaks into "
	"Rayleigh-Taylor fingers, and a quench lets it go as a fireball.",
	"Containment FFGL effect"                             // About
);

extern "C" const char* ContainmentBuildStamp()
{
	return "containment " CONTAINMENT_VERSION ", built " __DATE__ " " __TIME__;
}
