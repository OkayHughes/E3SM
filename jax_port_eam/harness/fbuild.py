"""Reusable f2py extraction-build helper (METHODOLOGY.md Tier-1).

Pattern: compile the real, unmodified E3SM Fortran sources plus
infrastructure stubs to objects with gfortran, then f2py-wrap a thin
hand-written driver (plain array interfaces — f2py handles those far
better than elemental/optional-argument routines) linking the objects.

Usage from a per-scheme build script:

    from fbuild import build_extension
    build_extension(
        name="eam_wv_sat_f",
        sources=[...ordered .F90 paths...],   # deps first
        driver="drivers/wv_sat_driver.F90",
        includes=[CAM_DIR],                   # for bfb_math.inc etc.
    )

Everything lands in jax_port_eam/build/<name>/ and the importable .so
is copied to jax_port_eam/fmod/ (kept out of git).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # jax_port_eam/
REPO = ROOT.parent                      # E3SM/
BUILD = ROOT / "build"
FMOD = ROOT / "fmod"

FFLAGS = ["-O2", "-fPIC", "-ffree-line-length-none",
          "-fallow-argument-mismatch", "-std=legacy"]

# canonical shared sources most schemes need, in dependency order
SHR_SOURCES = [
    REPO / "share/util/shr_kind_mod.F90",
    REPO / "share/util/shr_const_mod.F90",
    HERE / "stubs/infrastructure_stubs.F90",
    HERE / "stubs/physconst_stub.F90",
]


def build_extension(name, sources, driver, includes=(), fflags=()):
    """Compile `sources` (dependency order) to .o, then f2py-wrap
    `driver` linking them. Returns the path of the built extension."""
    bdir = BUILD / name
    bdir.mkdir(parents=True, exist_ok=True)
    FMOD.mkdir(exist_ok=True)

    flags = FFLAGS + list(fflags) + [f"-I{i}" for i in includes] \
        + ["-J", str(bdir), "-I", str(bdir)]
    objs = []
    for src in map(Path, sources):
        obj = bdir / (src.stem + ".o")
        cmd = ["gfortran", "-c", str(src), "-o", str(obj)] + flags
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        objs.append(str(obj))

    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep \
        + env.get("PATH", "")
    if "SDKROOT" not in env and sys.platform == "darwin":
        # macOS: gfortran needs the SDK to link
        sdk = subprocess.run(["xcrun", "--show-sdk-path"],
                             capture_output=True, text=True)
        if sdk.returncode == 0:
            env["SDKROOT"] = sdk.stdout.strip()
    cmd = [sys.executable, "-m", "numpy.f2py", "-c",
           str(Path(driver).resolve()), "-m", name,
           "--f2cmap", str(HERE / "f2py_f2cmap")] + objs \
        + [f"--f90flags={' '.join(FFLAGS + list(fflags))} -I{bdir} "
           + " ".join(f"-I{i}" for i in includes)]
    r = subprocess.run(cmd, cwd=bdir, capture_output=True, text=True,
                       env=env)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-4000:] + "\n" + r.stderr[-4000:])
        raise RuntimeError(f"f2py build of {name} failed")

    ext = next(bdir.glob(f"{name}*.so"))
    target = FMOD / ext.name
    # NEVER overwrite a signed dylib in place: macOS's kernel signature
    # cache goes stale and subsequent dlopen()s hang in uninterruptible
    # kernel waits. Copy to a temp name and rename (new inode, atomic).
    tmp = target.with_suffix(".so.new")
    shutil.copy2(ext, tmp)
    if target.exists():
        target.unlink()
    os.replace(tmp, target)
    print(f"built {target}")
    return target


def import_extension(name):
    """Import a previously built extension from fmod/."""
    if str(FMOD) not in sys.path:
        sys.path.insert(0, str(FMOD))
    return __import__(name)
