#!/usr/bin/env bash
#
# Every shader the plugin compiles, through a real GLSL compiler, before a host
# has to find out. Called by tools/verify.sh AND by CI, so the two cannot
# drift: a GitHub macOS runner cannot create an accelerated GL context, so in
# CI this is the only thing that looks at the GLSL at all. It is not a
# substitute for a driver -- Apple's Metal-backed GL can disagree with glslc --
# and that is checked on the dev Mac by every cttest GL check, and nowhere else.
#
#   tools/check-shaders.sh [path/to/cttest]
#
# The shaders are NOT scraped out of the source with a regex. The programs are
# assembled at run time from pieces in Shaders.cpp, so the text a regex could
# pull out is not the text that runs; `cttest --dump-shaders DIR` writes each
# program through the same SourceFor() the plugin compiles, and prints how many.
#
# Two checks:
#   - no GLSL reserved word is declared as an identifier (a driver may accept
#     `float sample;` and the next one refuses the whole program);
#   - every program compiles under glslc. --target-env=opengl4.5 with
#     -fauto-map-locations: glslc targets SPIR-V, which demands an explicit
#     layout( location ) on every uniform, a Vulkan rule rather than a GLSL one.
#     glslc is optional (`brew install shaderc`): without it that half skips.
#
set -uo pipefail
cd "$(dirname "$0")/.."

CTTEST="${1:-build/cttest}"
EXPECTED=14

if [ ! -x "$CTTEST" ]; then
	printf '   %s is not built\n' "$CTTEST"
	exit 1
fi

dir="$(mktemp -d)"
trap 'rm -rf "$dir"' EXIT

count=$("$CTTEST" --dump-shaders "$dir") || { printf '   cttest --dump-shaders failed\n'; exit 1; }
bad=0

# The count is asserted, not counted up to: a check that silently looks at
# fewer shaders than exist is worse than no check.
n=$(ls "$dir"/*.frag 2>/dev/null | wc -l | tr -d ' ')
if [ "$n" -ne "$EXPECTED" ] || [ "$count" != "$EXPECTED" ]; then
	printf '   %s programs written (cttest says %s), expected %d -- update EXPECTED with the list\n' "$n" "$count" "$EXPECTED"
	exit 1
fi

words="patch sample input output filter common active half layout flat smooth noperspective"
for word in $words; do
	if grep -nE "(float|int|uint|bool|vec[234]|ivec[234]|uvec[234]|mat[234]|Q|Side|Flux)[[:space:]]+$word[[:space:]]*[;=,)]" \
	            "$dir"/*.frag "$dir"/*.vert >/dev/null 2>&1; then
		printf '   "%s" is declared as an identifier and is a GLSL reserved word\n' "$word"
		bad=$(( bad + 1 ))
	fi
done
[ "$bad" -eq 0 ] && printf '   %d programs; no reserved word is used as an identifier\n' "$n"

if ! command -v glslc >/dev/null 2>&1; then
	printf '   glslc skipped: not installed (brew install shaderc)\n'
	exit "$bad"
fi
compiled=0
for shader in "$dir"/*.vert "$dir"/*.frag; do
	if glslc --target-env=opengl4.5 -fauto-map-locations "$shader" -o /dev/null 2>"$dir/err"; then
		compiled=$(( compiled + 1 ))
	else
		printf '   %s does not compile\n' "$(basename "$shader")"
		sed "s|$dir/||; s|^|      |" "$dir/err"
		bad=$(( bad + 1 ))
	fi
done
[ "$bad" -eq 0 ] && printf '   %d shaders (%d programs), all compile under glslc\n' "$compiled" "$n"
exit "$bad"
