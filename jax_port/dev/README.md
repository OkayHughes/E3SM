# Local EAMxx build environment (Docker)

Everything needed to build and test EAMxx (SCREAM) standalone on a laptop —
no supported HPC machine required. This is the environment the JAX-port test
harness ([../TEST_HARNESS_DESIGN.md](../TEST_HARNESS_DESIGN.md)) assumes for
Tiers 0–2.

The image bakes in the exact toolchain from
`components/eamxx/scripts/setup-copilot-env.sh` (Ubuntu 24.04, GNU compilers,
OpenMPI, NetCDF/PnetCDF, boost, yaml-cpp). The repo is bind-mounted, so the
source tree, the build tree (`components/eamxx/ctest-build/`), and the input
data (`../e3sm-inputdata`) all live on the host and survive container
deletion.

## One-time setup

```bash
# From the E3SM repo root, on the host
docker build -t scream-dev jax_port/dev

# Input data cache lives NEXT TO the repo (mounted at /work/e3sm-inputdata)
mkdir -p ../e3sm-inputdata

# Long-lived container; repo at /work/E3SM
docker run -d --name scream-dev \
    -v "$PWD":/work/E3SM \
    -v "$(cd .. && pwd)/e3sm-inputdata":/work/e3sm-inputdata \
    scream-dev tail -f /dev/null
```

Docker Desktop → Settings → Resources: give it **≥ 8 CPUs and ≥ 16 GB RAM**;
the debug C++ build is template-heavy and will thrash with less.

Submodules: the `jax_scream` working tree already has them checked out
(`git submodule update --init --recursive` on the host otherwise).

## Configure (once per build type)

```bash
docker exec -w /work/E3SM/components/eamxx scream-dev \
    ./scripts/test-all-eamxx -m copilot-testing -t dbg --config-only
```

This also downloads lookup tables / IC files from `web.lcrc.anl.gov` into
`/work/e3sm-inputdata` (needs network; one-time).

Build types: `dbg` → `ctest-build/copilot-testing/full_debug`,
`sp` → `full_sp_debug`, `opt` → `release`, `fpe` → `debug_nopack_fpe`.

## Build and test (the daily loop)

```bash
docker exec -w /work/E3SM/components/eamxx/ctest-build/copilot-testing/full_debug scream-dev \
    bash -c 'make -j$(nproc)'

docker exec -w /work/E3SM/components/eamxx/ctest-build/copilot-testing/full_debug scream-dev \
    ctest -R cld_fraction          # or: p3, shoc_standalone, --rerun-failed, -j$(nproc) ...
```

Interactive shell: `docker exec -it scream-dev bash`.
Container stopped/rebooted host: `docker start scream-dev`.

## Notes

- **Bind-mount I/O on macOS** is slower than container-local disk. Acceptable
  for incremental work; if full rebuilds hurt, move `ctest-build/` to a named
  volume (`-v scream-build:/work/E3SM/components/eamxx/ctest-build`).
- **Apple Silicon:** the image builds natively for linux/arm64 — no emulation.
  All apt packages and CPU JAX wheels exist for arm64.
- **Re-running `setup-copilot-env.sh` inside the container** is harmless
  (packages already present) but unnecessary. If it runs, its SSH probe fails
  in-container and it writes a `url.https://github.com/.insteadOf` rewrite
  into the *mounted repo's* `.git/config` — delete that line on the host if
  you use SSH remotes.
- **For the JAX harness** (later milestones): `pip install jax` in the
  container, and reconfigure with `EAMXX_ENABLE_PYTHON=ON` /
  `EAMXX_ENABLE_PYSCREAM=ON` (plus pybind11/nanobind, `pip install nanobind
  mpi4py`). Documented in TEST_HARNESS_DESIGN.md §5.
