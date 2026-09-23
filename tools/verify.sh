#!/usr/bin/env bash
#
# Everything, in the order that fails fastest.
#
# The build is universal on purpose. An arm64-only bundle builds and tests
# perfectly well here and then fails to load in an Intel Resolume, and the
# build log calls it a success either way -- so the architecture is checked
# with lipo, never with the log.
#
#     tools/verify.sh [BUILD_DIR]
#
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$REPO/build-verify}"

cd "$REPO"

step() { printf '\n\033[1m== %s\033[0m\n' "$1"; }
fail() { printf '\033[31mFAIL\033[0m %s\n' "$1"; exit 1; }

#---------------------------------------------------------------------------
# Every shader, exactly as the plugin compiles it.
#
# The programs are assembled at run time from pieces in Shaders.cpp, so the
# text a regex could pull out of the source is not the text that runs. The
# harness writes each program out through the same SourceFor() the plugin
# compiles (`cttest --dump-shaders DIR`) and both checks below read that.
#
# Needs a build, so it runs after one.
#---------------------------------------------------------------------------
reserved_words() {
	local dir="$1" words="patch sample input output filter common active half layout flat smooth noperspective"
	local bad=0 word
	for word in $words; do
		if grep -nE "(float|int|uint|bool|vec[234]|ivec[234]|uvec[234]|mat[234]|Q|Side|Flux)[[:space:]]+$word[[:space:]]*[;=,)]" \
		            "$dir"/*.frag "$dir"/*.vert >/dev/null 2>&1; then
			printf '   "%s" is declared as an identifier and is a GLSL reserved word\n' "$word"
			bad=$(( bad + 1 ))
		fi
	done
	[ "$bad" -eq 0 ] && printf '   none of the reserved words is used as an identifier\n'
	return "$bad"
}

# glslc targets SPIR-V; --target-env=opengl4.5 with -fauto-map-locations lets
# plain GLSL 4.10 through. Optional: `brew install shaderc`.
shaders_compile() {
	local dir="$1" bad=0 n=0 shader
	if ! command -v glslc >/dev/null 2>&1; then
		printf '   skipped: glslc not installed (brew install shaderc)\n'
		return 0
	fi
	for shader in "$dir"/*.vert "$dir"/*.frag; do
		[ -e "$shader" ] || continue
		n=$(( n + 1 ))
		if ! glslc --target-env=opengl4.5 -fauto-map-locations "$shader" -o /dev/null 2>"$dir/err"; then
			printf '   %s does not compile\n' "$( basename "$shader" )"
			sed "s|$dir/||; s|^|      |" "$dir/err"
			bad=$(( bad + 1 ))
		fi
	done
	if [ "$n" -eq 0 ]; then
		printf '   no shaders were dumped -- a check that looks at nothing is not a check\n'
		return 1
	fi
	[ "$bad" -eq 0 ] && printf '   %d shaders, all compile\n' "$n"
	return "$bad"
}

#---------------------------------------------------------------------------
step "Submodule"
#---------------------------------------------------------------------------
if [[ ! -f external/ffgl/CMakeLists.txt ]]; then
	fail "FFGL SDK missing -- run: git submodule update --init --recursive"
fi
echo "ok   FFGL SDK present at $(git -C external/ffgl rev-parse --short HEAD)"

#---------------------------------------------------------------------------
step "Build (universal)"
#---------------------------------------------------------------------------
cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$BUILD" -j"$(sysctl -n hw.ncpu)" >/dev/null
echo "ok   built"

shaders="$( mktemp -d )"
"$BUILD/cttest" --dump-shaders "$shaders" | sed 's/^/   /'

#---------------------------------------------------------------------------
step "GLSL reserved words"
#---------------------------------------------------------------------------
reserved_words "$shaders" || fail "a GLSL reserved word is used as an identifier"

#---------------------------------------------------------------------------
step "Shaders"
#---------------------------------------------------------------------------
shaders_compile "$shaders" || fail "a shader does not compile"
rm -rf "$shaders"

#---------------------------------------------------------------------------
step "Bundle"
#---------------------------------------------------------------------------
bundle="$BUILD/Containment.bundle"
binary="$bundle/Contents/MacOS/Containment"

[[ -f "$binary" ]] || fail "no binary at $binary"

# Universal. The failure this catches ships a plugin that simply does not
# appear in half the Resolume installs it is given to.
arches="$(lipo -archs "$binary")"
[[ "$arches" == *arm64* ]]  || fail "no arm64 slice (got: $arches)"
[[ "$arches" == *x86_64* ]] || fail "no x86_64 slice (got: $arches)"

# The entry point. A bundle whose registration got dropped by the linker still
# loads and still exports this -- the OBJECT-library note in CMakeLists.txt is
# what actually guards the registration; this catches a build that produced no
# module at all.
# Captured, then matched from a herestring -- never `nm ... | grep -q`.
# Under `set -o pipefail` a `grep -q` that finds its match exits immediately,
# the writer upstream takes SIGPIPE, and the PIPELINE reports failure even
# though the symbol is there. It is output-size dependent, so it fires on the
# bigger binary first and looks intermittent. A herestring is not a pipeline,
# so nothing can SIGPIPE.
symbols=$( nm -gU "$binary" 2>/dev/null || true )
grep -q '_plugMain' <<<"$symbols" || fail "plugMain not exported"

echo "ok   Containment: $arches, plugMain exported"

#---------------------------------------------------------------------------
step "Bundle metadata"
#---------------------------------------------------------------------------
plist="$bundle/Contents/Info.plist"
[[ -f "$plist" ]] || fail "no Info.plist in the bundle"

read_plist() { /usr/libexec/PlistBuddy -c "Print :$1" "$plist" 2>/dev/null || true; }

identifier="$( read_plist CFBundleIdentifier )"
executable="$( read_plist CFBundleExecutable )"
package="$( read_plist CFBundlePackageType )"
version="$( read_plist CFBundleVersion )"
declared="$( sed -n 's/^[[:space:]]*VERSION \([0-9.]*\)$/\1/p' CMakeLists.txt | head -1 )"

[[ "$identifier" == "com.stoatworks.ffgl.containment" ]] || fail "bundle id is '$identifier'"
[[ "$executable" == "Containment" ]] || fail "CFBundleExecutable is '$executable'"
[[ "$package" == "BNDL" ]] || fail "CFBundlePackageType is '$package', not BNDL"
# The version drifts across the manifest, the plist and the About header more
# often than anything else in the fleet, so all three are compared rather than
# any one of them trusted.
[[ "$version" == "$declared" ]] || fail "plist version '$version' != CMakeLists '$declared'"
grep -q "versionFallback = \"v$declared\"" source/StoatworksAbout.h \
	|| fail "StoatworksAbout.h's versionFallback is not v$declared"

echo "ok   $identifier, $package, v$version -- plist, CMakeLists and About agree"

#---------------------------------------------------------------------------
step "Code signature"
#---------------------------------------------------------------------------
# Ad hoc, which is what a local build gets. The release workflow signs and
# notarises properly; this only proves the bundle is well enough formed to be
# signed at all, which a malformed one is not.
codesign --force --sign - --timestamp=none "$bundle" >/dev/null 2>&1 \
	|| fail "the bundle could not be ad-hoc signed"
codesign --verify --deep --strict "$bundle" >/dev/null 2>&1 \
	|| fail "the ad-hoc signature does not verify"
echo "ok   ad-hoc signed and verified"

#---------------------------------------------------------------------------
step "Host view"
#---------------------------------------------------------------------------
# What a host actually reads out of the bundle: the id, the name and the type.
# The FFGL name field is char[ 16 ] and is NOT null-terminated, so a long name
# is truncated silently and nothing in this repo would ever notice -- only
# something that reads the bundle the way a host does.
#
# oxbow lives in the fleet, not here, so this is a skip rather than a failure
# when it is not to hand.
OXBOW="${OXBOW:-$HOME/Projects/resolume/oxbow/build/oxbow}"
if [[ -x "$OXBOW" ]]; then
	probe="$( "$OXBOW" probe "$bundle" 2>&1 || true )"
	printf '%s\n' "$probe" | sed 's/^/   /'
	grep -q 'CT01' <<<"$probe" || fail "oxbow did not read the id CT01"
	grep -q 'SW Containment' <<<"$probe" || fail "oxbow did not read the name 'SW Containment'"
	grep -qi 'effect' <<<"$probe" || fail "oxbow did not read the type as an effect"
	echo "ok   a host reads CT01 / SW Containment / effect"
else
	echo "   skipped: no oxbow at $OXBOW (set OXBOW=...)"
fi

#---------------------------------------------------------------------------
step "Checks"
#---------------------------------------------------------------------------
# Every claim the README makes, in the order the README makes them.
for check in briowu alfven conserve divb balance rt cusp frozen quench resist floors still glow state; do
	"$BUILD/cttest" --$check
done

#---------------------------------------------------------------------------
step "The harness drives the shipped shader"
#---------------------------------------------------------------------------
"$BUILD/cttest" --mutation

#---------------------------------------------------------------------------
step "Negative controls"
#---------------------------------------------------------------------------
# The checks above, run against a model that is deliberately wrong, and
# required to fail. A check that cannot fail is not a check.
"$BUILD/cttest" --negative

#---------------------------------------------------------------------------
step "Dead controls"
#---------------------------------------------------------------------------
# The only thing that catches a uniform whose name does not match the C++.
python3 tools/sweep.py --build "$(basename "$BUILD")"

#---------------------------------------------------------------------------
step "Cost"
#---------------------------------------------------------------------------
"$BUILD/cttest" --bench

printf '\n\033[32mall green\033[0m\n'
