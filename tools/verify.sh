#!/usr/bin/env bash
#
# Everything that can be checked without a host, in one go.
#
#     tools/verify.sh [BUILD_DIR]        (default build-verify)
#
# Each step answers a question none of the others can:
#
#   submodule     the FFGL SDK is there, at the fleet's pin (b1afaf9).
#   build         a FRESH universal Release build. Not the dev build: CMake
#                 latches the architecture list at the first target, so the
#                 only build worth measuring is one configured from nothing.
#   shaders       every program the plugin compiles, as it compiles it: no
#                 reserved word as an identifier, and glslc accepts it.
#   offline       the checks that need no GL -- names, presets, the Brio-Wu
#                 reference, the coils' vacuum field -- and their negative
#                 controls. This is what CI runs.
#   physics       every GL check, TWICE: at each check's own raster, and with
#                 every rig at 320x180, CI's. The grid is the check's either way
#                 (Detail cells on the short side, the aspect the check asked
#                 for); the second pass moves only the light's path to the
#                 raster. A check that held at one raster was fitted to it.
#   mutation      one character of the shipped limiter changed: checks fail.
#   negative      every check against a deliberately wrong model: each fails.
#   pipe          the fleet's --pipe frame format: a partial frame at EOF ends
#                 the stream, a cue naming no control is refused, and a reader
#                 that hangs up ends the run with exit 1, not SIGPIPE's 141.
#   sweep         does every control change the picture.
#   bundle        lipo (universal), plugMain, the plist, an ad-hoc signature.
#   oxbow         a real FFGL host loads the bundle and reports the name, id
#                 and type it sees -- the name field is not null-terminated and
#                 a host truncates silently past 16 characters.
#   bench         the render cost, for the record. Not pass/fail.
#
# Nothing stops at the first failure: every step runs, and the summary at the
# end counts what failed.
#
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$REPO/build-verify}"
cd "$REPO"

failures=0
passes=0
step() { printf '\n\033[1m== %s\033[0m\n' "$1"; }
pass() { printf '   \033[32mok\033[0m   %s\n' "$1"; passes=$(( passes + 1 )); }
fail() { printf '   \033[31mFAIL\033[0m %s\n' "$1"; failures=$(( failures + 1 )); }

#---------------------------------------------------------------------------
step "submodule"
#---------------------------------------------------------------------------
if [[ ! -f external/ffgl/CMakeLists.txt ]]; then
	fail "FFGL SDK missing -- run: git submodule update --init --recursive"
	printf '\n%d passed, %d FAILED\n' "$passes" "$failures"
	exit 1
fi
pin="$( git -C external/ffgl rev-parse --short=7 HEAD )"
if [[ "$pin" == "b1afaf9" ]]; then
	pass "FFGL SDK at $pin, the fleet's pin"
else
	fail "FFGL SDK at $pin, not the fleet's b1afaf9"
fi

#---------------------------------------------------------------------------
step "build (fresh universal Release, $(basename "$BUILD"))"
#---------------------------------------------------------------------------
rm -rf "$BUILD"
buildlog="$( mktemp )"
if cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release >"$buildlog" 2>&1 \
	&& cmake --build "$BUILD" -j"$(sysctl -n hw.ncpu)" >>"$buildlog" 2>&1; then
	pass "built ($( grep -c 'warning:' "$buildlog" ) warnings, all in the vendored SDK: $( grep 'warning:' "$buildlog" | grep -vc 'external/ffgl' ) outside it)"
	rm -f "$buildlog"
else
	tail -30 "$buildlog" | sed 's/^/      /'
	fail "the build failed"
	printf '\n%d passed, %d FAILED\n' "$passes" "$failures"
	exit 1
fi
CTTEST="$BUILD/cttest"

#---------------------------------------------------------------------------
step "shaders"
#---------------------------------------------------------------------------
if out=$( tools/check-shaders.sh "$CTTEST" ); then
	printf '%s\n' "$out"
	pass "every program: no reserved word, glslc compiles it"
else
	printf '%s\n' "$out"
	fail "a shader is not clean"
fi

#---------------------------------------------------------------------------
step "offline"
#---------------------------------------------------------------------------
if out=$( "$CTTEST" --offline 2>&1 ); then
	pass "cttest --offline: names, presets, reference, vacuum; $( grep -o 'negative controls: .*' <<<"$out" )"
else
	printf '%s\n' "$out" | grep -E 'FAIL' | sed 's/^/      /'
	fail "cttest --offline"
fi

#---------------------------------------------------------------------------
step "physics, at each check's own raster and at 320x180"
#---------------------------------------------------------------------------
# Every claim the README makes, in the order it makes them.
checks="briowu alfven conserve divb balance rt cusp frozen quench resist floors still glow state open"
for size in own 320x180; do
	for check in $checks; do
		args=( "--$check" )
		[[ "$size" != own ]] && args+=( --size "$size" )
		if out=$( "$CTTEST" "${args[@]}" 2>&1 ); then
			pass "--$check @ $size: $( grep -c '^  ok' <<<"$out" ) assertions"
		else
			fail "--$check @ $size"
			printf '%s\n' "$out" | grep -E 'FAIL ' | sed 's/^/      /'
		fi
	done
done

#---------------------------------------------------------------------------
step "mutation: the harness drives the shipped shader"
#---------------------------------------------------------------------------
if out=$( "$CTTEST" --mutation 2>&1 ); then
	pass "$( grep -o 'against the mutated shader.*' <<<"$out" )"
else
	fail "cttest --mutation"
	printf '%s\n' "$out" | tail -3 | sed 's/^/      /'
fi

#---------------------------------------------------------------------------
step "negative controls"
#---------------------------------------------------------------------------
# The checks above, run against a model that is deliberately wrong, and
# required to fail. A check that cannot fail is not a check.
if out=$( "$CTTEST" --negative 2>&1 ); then
	pass "$( grep -o 'negative controls: .*' <<<"$out" )"
else
	fail "a negative control passed against a wrong model"
	printf '%s\n' "$out" | grep -E 'PASSED against' | sed 's/^/      /'
fi

#---------------------------------------------------------------------------
step "pipe"
#---------------------------------------------------------------------------
# Two and a half frames in must be exactly two frames out and a clean exit --
# a partial frame is the end of the stream, never a frame.
frame=$(( 64 * 36 * 4 ))
raw=$( mktemp ); cues=$( mktemp )
head -c $(( frame * 5 / 2 )) /dev/zero > "$raw"
got=$( "$CTTEST" --pipe --size 64x36 < "$raw" 2>/dev/null | wc -c | tr -d ' ' )
status=${PIPESTATUS[0]}
if [ "$status" -eq 0 ] && [ "$got" = "$(( frame * 2 ))" ]; then
	pass "2.5 frames in, exactly 2 frames out, clean exit"
else
	fail "2.5 frames in gave $got bytes out (want $(( frame * 2 ))), exit $status"
fi
# A cue naming no parameter must be refused rather than silently doing nothing
# to a take. Read from a file, not a pipe: a writer killed by SIGPIPE would fail
# the pipeline whatever cttest did, and the refusal would pass for the wrong
# reason.
printf '0 No Such Control 0.5\n' > "$cues"
"$CTTEST" --pipe --size 64x36 --script "$cues" < "$raw" >/dev/null 2>&1
status=$?
if [ "$status" -eq 2 ]; then
	pass "a cue naming no parameter is refused (exit 2)"
else
	fail "a cue naming no parameter gave exit $status, not 2"
fi
printf '0 Field 0.25\n1 Field 0.5\n0 Ignite 0\n1 Ignite 1\n' > "$cues"
got=$( "$CTTEST" --pipe --size 64x36 --script "$cues" < "$raw" 2>/dev/null | wc -c | tr -d ' ' )
if [ "$got" = "$(( frame * 2 ))" ]; then
	pass "a cue sheet of real controls is accepted and still gives 2 frames"
else
	fail "a cue sheet of real controls gave $got bytes"
fi
# The picture goes through: a mid-grey clip with Mix 0 comes back mid-grey.
python3 -c "import sys; sys.stdout.buffer.write(bytes([128,128,128,255]) * (64*36*2))" > "$raw"
if python3 - "$CTTEST" "$raw" <<'PY'
import subprocess, sys
out = subprocess.run([sys.argv[1], "--pipe", "--size", "64x36", "--set", "Mix=0"],
                     stdin=open(sys.argv[2], "rb"), capture_output=True).stdout
sys.exit(0 if out == open(sys.argv[2], "rb").read() else 1)
PY
then
	pass "with Mix 0 the frames come back byte for byte"
else
	fail "with Mix 0 the piped frames are not the input"
fi
# A reader that hangs up early (`| head -c 1`, ffmpeg dying) must end the run
# with exit 1 and a message, not SIGPIPE's silent 141.
head -c $(( frame * 20 )) /dev/zero > "$raw"
"$CTTEST" --pipe --size 64x36 < "$raw" 2>/dev/null | head -c 1 >/dev/null
status=${PIPESTATUS[0]}
if [ "$status" -eq 1 ]; then
	pass "a closed stdout ends the run with exit 1, not SIGPIPE"
else
	fail "a closed stdout gave exit $status, not 1"
fi
rm -f "$raw" "$cues"

#---------------------------------------------------------------------------
step "sweep"
#---------------------------------------------------------------------------
# The only thing that catches a uniform whose name does not match the C++.
if out=$( python3 tools/sweep.py --binary "$CTTEST" 2>&1 ); then
	pass "$( printf '%s\n' "$out" | tail -1 )"
else
	printf '%s\n' "$out" | sed 's/^/      /'
	fail "a control is dead"
fi

#---------------------------------------------------------------------------
step "bundle"
#---------------------------------------------------------------------------
bundle="$BUILD/Containment.bundle"
binary="$bundle/Contents/MacOS/Containment"
if [[ -f "$binary" ]]; then
	# Universal. The failure this catches ships a plugin that simply does not
	# appear in half the Resolume installs it is given to.
	arches="$( lipo -archs "$binary" )"
	if [[ "$arches" == *arm64* && "$arches" == *x86_64* ]]; then
		pass "universal: $arches"
	else
		fail "not universal: $arches"
	fi
	# Captured, then matched from a herestring -- never `nm ... | grep -q`.
	# Under `set -o pipefail` a `grep -q` that finds its match exits at once,
	# nm takes SIGPIPE, and the PIPELINE reports failure though the symbol is
	# there.
	symbols=$( nm -gU "$binary" 2>/dev/null || true )
	if grep -q '_plugMain' <<<"$symbols"; then pass "plugMain exported"; else fail "plugMain not exported"; fi

	plist="$bundle/Contents/Info.plist"
	read_plist() { /usr/libexec/PlistBuddy -c "Print :$1" "$plist" 2>/dev/null || true; }
	identifier="$( read_plist CFBundleIdentifier )"
	executable="$( read_plist CFBundleExecutable )"
	package="$( read_plist CFBundlePackageType )"
	version="$( read_plist CFBundleVersion )"
	declared="$( sed -n 's/^[[:space:]]*VERSION \([0-9.]*\)$/\1/p' CMakeLists.txt | head -1 )"
	# The version drifts across the manifest, the plist and the About header
	# more often than anything else in the fleet, so all three are compared.
	if [[ "$identifier" == "com.stoatworks.ffgl.containment" && "$executable" == "Containment" && "$package" == "BNDL" \
		&& "$version" == "$declared" ]] && grep -q "versionFallback = \"v$declared\"" source/StoatworksAbout.h; then
		pass "$identifier, $package, v$version -- plist, CMakeLists and About agree"
	else
		fail "plist: id '$identifier', executable '$executable', type '$package', version '$version' (CMakeLists $declared)"
	fi
	# Ad hoc, which is what a local build gets; the release workflow signs and
	# notarises properly. This proves the bundle is well enough formed to sign.
	if codesign --force --sign - --timestamp=none "$bundle" >/dev/null 2>&1 \
		&& codesign --verify --deep --strict "$bundle" >/dev/null 2>&1; then
		pass "ad-hoc signed and verified"
	else
		fail "the bundle does not take an ad-hoc signature"
	fi
else
	fail "no binary at $binary"
fi

#---------------------------------------------------------------------------
step "oxbow"
#---------------------------------------------------------------------------
# What a host actually reads out of the bundle. oxbow lives in the fleet, not
# here, so it is a skip rather than a failure when it is not to hand.
OXBOW="${OXBOW:-$HOME/Projects/resolume/oxbow/build/oxbow}"
if [[ -x "$OXBOW" ]]; then
	probe="$( "$OXBOW" probe "$bundle" 2>&1 || true )"
	printf '%s\n' "$probe" | sed 's/^/      /'
	if grep -q 'CT01' <<<"$probe" && grep -q 'SW Containment' <<<"$probe" && grep -qi 'effect' <<<"$probe"; then
		pass "a host reads CT01 / SW Containment / effect"
	else
		fail "oxbow did not read CT01 / SW Containment / effect"
	fi
else
	printf '   skipped: no oxbow at %s (set OXBOW=...)\n' "$OXBOW"
fi

#---------------------------------------------------------------------------
step "bench (for the record)"
#---------------------------------------------------------------------------
"$CTTEST" --bench | sed 's/^/   /'

printf '\n'
if [ "$failures" -eq 0 ]; then
	printf '\033[32mall green\033[0m: %d passed, 0 failed\n' "$passes"
	exit 0
fi
printf '\033[31m%d FAILED\033[0m, %d passed\n' "$failures" "$passes"
exit 1
