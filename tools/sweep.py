#!/usr/bin/env python3
"""Move every parameter and fail if any of them made no difference to the frame.

**This is the only thing in the repo that catches a dead control**, and it is
not a theoretical risk. A GLSL uniform whose name does not match the C++ is
silently ignored -- `glGetUniformLocation` returns -1 and `glUniform` on -1 is
a documented no-op -- so a slider can be stone dead while everything compiles,
links, loads and renders. The measurement checks will not catch it either: they
only ever exercise the settings they set themselves.

## The context table

Most of Containment's controls act on the default look straight away: the
ball is ignited on the first frame, the coils are turning, the drive stirs.
A few need something else to be true first, supplied here as raw harness
arguments and settings and nothing more:

- the two audio amounts need audio (`--beat` feeds a beat every half second);
- Field Line Count needs Field Lines up; Ramp needs Temperature Tint up.

An entry whose absence changes nothing is worse than no entry: it makes the
sweep look more careful than it is. Each one here was checked by removing it.

Usage::

    tools/sweep.py [--build BUILD_DIR] [--verbose]
"""

import argparse
import pathlib
import subprocess
import sys
import tempfile
import zlib

REPO = pathlib.Path(__file__).resolve().parent.parent

# Applied to every render: nothing -- the defaults are a live plasma.
BASE = []

# What else has to be true for a parameter to be able to do anything.
BEAT = ["--beat"]
CONTEXT = {
    "Audio Heat": BEAT,
    "Audio Pellets": BEAT,
}
SETTINGS = {
    "Field Line Count": ["Field Lines=0.8"],
    "Ramp": ["Temperature Tint=1"],
}

# The values every non-option parameter is swept across. The awkward numbers
# are load-bearing: a parameter the picture is periodic in can land 0, 0.5 and
# 1 on pixel-identical frames and report a working slider as dead. 0.137 and
# 0.611 are not rational multiples of anything swept here.
SWEEP_VALUES = [0.0, 0.137, 0.611, 1.0]

# Option parameters are swept across their elements instead. --list reports a
# parameter's kind but not its element count, so these track the enums in
# Controls.h by hand.
OPTION_RANGE = {
    "Profile": 2,
    "Poles": 5,
    "Quench": 2,
    "Boundary": 2,
    "Detail": 4,
    "Ramp": 2,
    "View": 7,
}

# Parameter kinds with no scalar worth sweeping.
SKIP_KINDS = {"buffer", "event", "text"}


def read_png(path):
    """Decode a PNG to raw bytes. Enough of the format for our own writer's
    output -- 8-bit RGBA, one IDAT, filter 0 on every row."""
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")

    pos = 8
    idat = b""
    while pos < len(data):
        length = int.from_bytes(data[pos:pos + 4], "big")
        kind = data[pos + 4:pos + 8]
        body = data[pos + 8:pos + 8 + length]
        if kind == b"IDAT":
            idat += body
        pos += 12 + length

    return zlib.decompress(idat)


def render(harness, out, settings, raw, verbose):
    args = [str(harness), "--out", str(out), "--size", "480x270", "--frames", "90"]
    args += raw
    for setting in settings:
        args += ["--set", setting]

    if verbose:
        print("   ", " ".join(args))

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"cttest failed: {result.stderr.strip()}")

    return read_png(out)


def parameters(harness):
    """Name and kind of every parameter, in declaration order."""
    result = subprocess.run([str(harness), "--list"], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"cttest --list failed: {result.stderr.strip()}")

    found = []
    for line in result.stdout.splitlines()[1:]:
        # id, name (may contain spaces), kind, default
        parts = line.split()
        if len(parts) < 4:
            continue
        kind = parts[-2]
        name = " ".join(parts[1:-2])
        found.append((name, kind))

    # The About block -- a text line then one browser button per link -- is
    # declared last and never touches a pixel, so sweeping it only buries a
    # real dead control. Truncating at "About" rather than naming each button
    # keeps this right as links are added; publishing a user guide adds one.
    for index, (entry_name, _kind) in enumerate(found):
        if entry_name == "About":
            return found[:index]

    return found


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", default="build")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    harness = REPO / args.build / "cttest"
    if not harness.exists():
        print(f"no cttest at {harness} -- build first", file=sys.stderr)
        return 2

    dead = []
    checked = 0

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)

        for name, kind in parameters(harness):
            if kind in SKIP_KINDS:
                continue

            extra = CONTEXT.get(name, [])
            base = BASE + SETTINGS.get(name, [])

            if name in OPTION_RANGE:
                values = [float(i) for i in range(OPTION_RANGE[name])]
            else:
                values = SWEEP_VALUES

            frames = []
            for value in values:
                out = tmp / "sweep.png"
                frames.append(
                    render(harness, out, base + [f"{name}={value}"], extra, args.verbose))

            checked += 1
            if all(f == frames[0] for f in frames[1:]):
                dead.append(f"{name} ({kind})")
                print(f"  DEAD {name}")
            elif args.verbose:
                print(f"  ok   {name}")

    print()
    if dead:
        print(f"sweep: {checked} parameters, {len(dead)} made no difference:")
        for entry in dead:
            print(f"  - {entry}")
        return 1

    print(f"sweep: {checked} parameters, all live")
    return 0


if __name__ == "__main__":
    sys.exit(main())
