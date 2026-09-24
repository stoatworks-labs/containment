"""The demo's shaders must be the plugin's shaders, character for character.

    python3 demo/tools/check_shaders.py

Called from `tools/verify.sh`. Exit code 1 means the two copies have drifted.

------------------------------------------------------------------- why

`demo/plugin.js` holds a copy of every piece of GLSL in `source/Shaders.cpp`.
Two copies drift -- quietly, because a plasma that looks *plausible* looks
exactly like the right one. The page's whole claim is that it runs the
plugin's own solver, so something has to enforce it. Nothing else can:
`tools/check-shaders.sh` hands the programs `cttest --dump-shaders` writes to
glslc and never looks at the page, and `cttest` drives the real plugin class
and has no idea this page exists.

------------------------------------------------------------------- what it does

Two things.

1. **The pieces.** Every program is assembled at run time by `SourceFor()`,
   which joins raw-string pieces (`kVersion`, `kCommon`, `kRiemann`,
   `kForces`, one `...Main` per pass). The page carries the same pieces as
   backtick literals. Each `R"( ... )"` body is pulled out of the C++ and the
   matching literal out of `plugin.js`, and compared exactly -- no whitespace
   normalisation, no comment stripping: a comment updated on one side and not
   the other is drift worth catching, because the comments carry the
   reasoning. `kVersion` is an ordinary C string, so its `\\n` is decoded.

2. **The assembly.** Each `case Program::X: return { "name", vertex, join( {
   ... } ) }` in `SourceFor()` must be matched by a `name: [ ... ]` row of
   `PROGRAM_PIECES` in plugin.js naming the same pieces in the same order,
   and the vertex shader must be `kVersion + kVertex` on both sides. The
   right pieces joined in the wrong order would be a different shader.

The one transformation on the JS side is a decode, not a normalisation. The
emission main quotes a command in a comment with backticks -- `cttest --divb`
-- and a backtick cannot appear raw in a JavaScript template literal, so
plugin.js escapes it as \\`. This unescapes that and *rejects any other
backslash and any `${`*; neither occurs in the C++ raw strings, so either
could only be somebody hiding a difference.

------------------------------------------------------------------- what it cannot

Nothing here checks the *ported* half: the conversions of `Controls.cpp`, the
CPU half of `Physics.cpp` (`MakeCoils`, `VacuumPotential`, `ChooseGrid`,
`Pcg`, `Drive`), `Presets.h` and `P()`, and the frame sequence of
`ProcessOpenGL` -- the substep plan, the events, `EnsureBuffers`, `Glow`'s
weights, `Potential`'s V-cycle. Those are a hand translation in plugin.js, and
only a reader can tell whether they still agree. When you change one of
those, change it there too.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

# JS constant, C++ symbol, and whether the C++ is a raw string.
SHADERS = [
    ("VERSION", "kVersion", False),
    ("VERTEX", "kVertex", True),
    ("COMMON", "kCommon", True),
    ("RIEMANN", "kRiemann", True),
    ("FORCES", "kForces", True),
    ("REDUCE_MAIN", "kReduceMain", True),
    ("CLOCK_MAIN", "kClockMain", True),
    ("PREDICT_MAIN", "kPredictMain", True),
    ("FLUX_MAIN", "kFluxMain", True),
    ("UPDATE_MAIN", "kUpdateMain", True),
    ("IGNITE_MAIN", "kIgniteMain", True),
    ("SOURCES_MAIN", "kSourcesMain", True),
    ("EMISSION_MAIN", "kEmissionMain", True),
    ("DOWNSAMPLE_MAIN", "kDownsampleMain", True),
    ("BLUR_MAIN", "kBlurMain", True),
    ("POISSON_COMMON", "kPoissonCommon", True),
    ("RESIDUAL_MAIN", "kResidualMain", True),
    ("SMOOTH_MAIN", "kSmoothMain", True),
    ("PROLONG_MAIN", "kProlongMain", True),
    ("COMPOSITE_MAIN", "kCompositeMain", True),
]
JS_NAME = {symbol: name for name, symbol, _ in SHADERS}

# The C string escapes kVersion could plausibly use. Anything else is refused.
C_ESCAPES = {"n": "\n", "t": "\t", "\\": "\\", '"': '"'}


def from_cpp(source, symbol, raw):
    if raw:
        match = re.search(r'const char\* const ' + symbol + r' = R"\((.*?)\)";', source, re.S)
        return None if match is None else match.group(1)
    match = re.search(r'const char\* const ' + symbol + r' = "((?:[^"\\]|\\.)*)";', source)
    if match is None:
        return None
    out, text, i = [], match.group(1), 0
    while i < len(text):
        if text[i] == "\\":
            if text[i + 1] not in C_ESCAPES:
                return None
            out.append(C_ESCAPES[text[i + 1]])
            i += 2
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def from_js(source, name):
    match = re.search(r'^const ' + name + r' = `(.*?)`;$', source, re.S | re.M)
    if match is None:
        return None, None
    body = match.group(1)
    stray = re.search(r"\\(?!`)", body)
    if stray is not None:
        upto = body[: stray.start()]
        return None, f"backslash that is not an escaped backtick, at line {upto.count(chr(10)) + 1}"
    if "${" in body:
        upto = body[: body.index("${")]
        return None, f"template interpolation, at line {upto.count(chr(10)) + 1}"
    return body.replace("\\`", "`"), None


def pieces_list(text):
    return [p.strip() for p in text.split(",") if p.strip()]


def main():
    with open(os.path.join(REPO, "source", "Shaders.cpp")) as handle:
        cpp = handle.read()
    with open(os.path.join(REPO, "demo", "plugin.js")) as handle:
        js = handle.read()

    problems = 0
    for name, symbol, raw in SHADERS:
        cpp_text = from_cpp(cpp, symbol, raw)
        js_text, complaint = from_js(js, name)
        if cpp_text is None:
            print(f"FAIL  {symbol} not found in source/Shaders.cpp")
            problems += 1
            continue
        if complaint is not None:
            print(f"FAIL  {name} in demo/plugin.js has a {complaint}")
            problems += 1
            continue
        if js_text is None:
            print(f"FAIL  {name} not found in demo/plugin.js")
            problems += 1
            continue
        if cpp_text == js_text:
            print(f"ok    {name:<16} matches {symbol} ({len(cpp_text)} chars)")
            continue
        problems += 1
        print(f"FAIL  {name} has drifted from {symbol}")
        a_lines, b_lines = cpp_text.splitlines(), js_text.splitlines()
        for i in range(max(len(a_lines), len(b_lines))):
            a = a_lines[i] if i < len(a_lines) else "<missing>"
            b = b_lines[i] if i < len(b_lines) else "<missing>"
            if a != b:
                print(f"        first difference at line {i + 1}")
                print(f"          C++: {a}")
                print(f"          js : {b}")
                break
        else:
            print("        the lines agree; the difference is a trailing newline")

    # Every piece the C++ declares must be carried and compared.
    declared = set(re.findall(r"const char\* const (k\w+) =", cpp))
    listed = {symbol for _, symbol, _ in SHADERS}
    for symbol in sorted(declared - listed):
        print(f"FAIL  {symbol} is in source/Shaders.cpp and this check does not compare it")
        problems += 1

    # The assembly: SourceFor() against PROGRAM_PIECES.
    if re.search(r"const std::string vertex = join\( \{ kVersion, kVertex \} \);", cpp) is None:
        print("FAIL  SourceFor() no longer builds the vertex shader as kVersion + kVertex")
        problems += 1
    elif "const VERTEX_SOURCE = VERSION + VERTEX;" not in js:
        print("FAIL  demo/plugin.js does not build VERTEX_SOURCE as VERSION + VERTEX")
        problems += 1
    else:
        print("ok    vertex           assembled as kVersion + kVertex")

    cases = re.findall(r'case Program::\w+:\s*return \{ "(\w+)", vertex, join\( \{ ([^}]*) \} \) \};', cpp)
    table = re.search(r"^const PROGRAM_PIECES = \{\n(.*?)^\};", js, re.S | re.M)
    rows = {} if table is None else dict(re.findall(r"^  (\w+): \[([^\]]*)\],$", table.group(1), re.M))
    if not cases:
        print("FAIL  no `case Program::...: return { ... join( { ... } ) }` found in SourceFor()")
        problems += 1
    for name, pieces in cases:
        want = [JS_NAME.get(p, f"<{p}?>") for p in pieces_list(pieces)]
        have = pieces_list(rows[name]) if name in rows else None
        if have is None:
            print(f"FAIL  program {name} is not in PROGRAM_PIECES")
            problems += 1
        elif have != want:
            print(f"FAIL  program {name}: plugin.js joins {have}, SourceFor() joins {want}")
            problems += 1
        else:
            print(f"ok    {name:<16} assembled as SourceFor() does ({len(want)} pieces)")
    for name in sorted(set(rows) - {n for n, _ in cases}):
        print(f"FAIL  PROGRAM_PIECES has {name}, which SourceFor() does not build")
        problems += 1

    print()
    if problems:
        print(f"{problems} problem(s) -- copy the C++ across, do not edit plugin.js by hand")
        return 1
    print(f"all {len(SHADERS)} shader pieces and {len(cases)} programs are identical to the plugin's")
    return 0


if __name__ == "__main__":
    sys.exit(main())
