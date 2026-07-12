# Development container for building/testing EAMxx (SCREAM) standalone,
# used by the JAX-port workflow (see jax_port/TEST_HARNESS_DESIGN.md).
#
# Mirrors components/eamxx/scripts/setup-copilot-env.sh: same apt package set,
# same multiarch symlink fix, same Python packages — baked into the image so
# containers are disposable. The repo is bind-mounted at /work/E3SM (see
# jax_port/dev/README.md), so source, build tree, and input data live on the
# host and survive container deletion.

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Package list from setup-copilot-env.sh (install_with_apt), plus sudo and
# openssh-client which the script itself expects to find.
RUN apt-get update && apt-get install -y --no-install-recommends \
        sudo git ca-certificates curl wget openssh-client libxml2-utils \
        gfortran gcc g++ cmake make m4 \
        python3 python3-pip python3-venv python3-dev \
        libopenmpi-dev openmpi-bin \
        libnetcdf-dev libnetcdff-dev \
        libboost-dev libyaml-cpp-dev \
        libblas-dev liblapack-dev \
        libcurl4-openssl-dev zlib1g-dev \
        libhdf5-openmpi-dev \
        perl \
    && rm -rf /var/lib/apt/lists/*

# Multiarch fix from setup-copilot-env.sh: Scorpio's Find modules only search
# /usr/lib with NO_DEFAULT_PATH, but Debian/Ubuntu put libs in lib/<triplet>/.
# blas/lapack added to the script's list: copilot-testing.cmake's fallback
# hardcodes the x86_64 multiarch dir, so on arm64 we route it through
# BLAS_ROOT=/usr (below), which expects /usr/lib/libblas.so.
RUN arch_dir="/usr/lib/$(gcc -dumpmachine)" && \
    for lib in libnetcdf.so libnetcdff.so libhdf5.so libhdf5_hl.so \
               libblas.so liblapack.so; do \
        if [ -f "$arch_dir/$lib" ] && [ ! -f "/usr/lib/$lib" ]; then \
            ln -sf "$arch_dir/$lib" "/usr/lib/$lib"; \
        fi; \
    done

# PnetCDF from source: Ubuntu 24.04's libpnetcdf-dev ships its Fortran .mod
# for a newer gfortran than the distro default, so Scorpio's Fortran interface
# can't use it (and E3SM's Scorpio requires PnetCDF and/or NetCDF-Fortran to
# be fully usable). Building against the container's exact mpif90 fixes this.
RUN curl -fsSL https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz \
        | tar -xz -C /tmp && \
    cd /tmp/pnetcdf-1.12.3 && \
    ./configure --prefix=/usr/local --enable-shared MPICC=mpicc MPICXX=mpicxx MPIF77=mpif77 MPIF90=mpif90 \
        >/tmp/pnetcdf-configure.log && \
    make -j"$(nproc)" >/tmp/pnetcdf-make.log && \
    make install >/tmp/pnetcdf-install.log && \
    ldconfig && \
    rm -rf /tmp/pnetcdf-1.12.3

# E3SM's FindNetCDF.cmake reads NETCDF_PATH/PNETCDF_PATH; Scorpio's own
# LibFind-based modules read NETCDF/NETCDF_FORTRAN/PNETCDF. /usr works for
# NetCDF because nc-config/nf-config live in /usr/bin and netcdf.mod in
# /usr/include; PnetCDF points at the /usr/local build above.
ENV NETCDF_PATH=/usr \
    NETCDF=/usr \
    NETCDF_FORTRAN=/usr \
    PNETCDF_PATH=/usr/local \
    PNETCDF=/usr/local \
    BLAS_ROOT=/usr

# Ubuntu 24.04 Python is PEP-668 "externally managed"; this container is
# single-purpose, so let pip install into the system environment. The env var
# also makes the pip step inside setup-copilot-env.sh work if it is re-run.
ENV PIP_BREAK_SYSTEM_PACKAGES=1
# Base EAMxx deps plus the embedded-Python/JAX harness deps: pybind11 for
# EAMXX_ENABLE_PYTHON, jax/pytest for the scream_jax swap tests (Tier 2 in
# jax_port/TEST_HARNESS_DESIGN.md).
RUN python3 -m pip install --no-cache-dir \
        psutil pyyaml netCDF4 packaging \
        "pybind11[global]" jax pytest

# OpenMPI quality-of-life inside containers: ctest launches mpiexec as root,
# and the CMA/vader single-copy path needs ptrace caps Docker doesn't grant.
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# The bind-mounted repo is owned by a host UID; without this git refuses to
# operate on it ("dubious ownership").
RUN git config --system --add safe.directory '*'

# Input data downloaded by configure lands here (on the bind mount, so it is
# downloaded once and shared across containers).
ENV SCREAM_INPUT_ROOT=/work/e3sm-inputdata

WORKDIR /work/E3SM
CMD ["/bin/bash"]
