#!/usr/bin/env bash
# Builds SLATE (with its BLAS++ and LAPACK++ submodules) into
# third_party/install, for use with `bazel build --define slate=1`.
#
# Needs: CMake, an MPI compiler (OpenMPI), OpenBLAS, and for the GPU backend
# the CUDA toolkit. Usage:
#   tools/build_slate.sh            # CUDA backend if nvcc is found
#   GPU=none tools/build_slate.sh   # CPU only
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TAG="${SLATE_TAG:-v2024.10.29}"
SRC="$ROOT/third_party/src/slate"
BUILD="$ROOT/third_party/build/slate"
PREFIX="$ROOT/third_party/install"
CUDA_HOME="${CUDA_HOME:-/usr/local/cuda-12.6}"
JOBS="${JOBS:-$(nproc)}"

if [[ "${GPU:-auto}" == "auto" ]]; then
  if [[ -x "$CUDA_HOME/bin/nvcc" ]]; then GPU=cuda; else GPU=none; fi
fi

command -v mpicxx >/dev/null || { echo "mpicxx not found (install libopenmpi-dev)"; exit 1; }

if [[ ! -d "$SRC" ]]; then
  git clone --depth 1 --branch "$TAG" --recurse-submodules --shallow-submodules \
    https://github.com/icl-utk-edu/slate.git "$SRC"
fi

CMAKE_ARGS=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
  -DCMAKE_CXX_COMPILER=mpicxx
  -DBUILD_SHARED_LIBS=ON
  -Dblas=openblas
  -Dbuild_tests=no
  -Dc_api=no
  -Dfortran_api=no
)
if [[ "$GPU" == "cuda" ]]; then
  export PATH="$CUDA_HOME/bin:$PATH"
  CMAKE_ARGS+=(-Dgpu_backend=cuda -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCH:-86}"
               -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc")
else
  CMAKE_ARGS+=(-Dgpu_backend=none)
fi

cmake -S "$SRC" -B "$BUILD" "${CMAKE_ARGS[@]}"
cmake --build "$BUILD" -j "$JOBS"
cmake --install "$BUILD"

# Record the MPI and CUDA flags for third_party/slate/slate.BUILD.
{
  echo "MPI_INCLUDE=$(mpicxx --showme:incdirs | awk '{print $1}')"
  echo "GPU=$GPU"
  echo "CUDA_HOME=$CUDA_HOME"
} > "$PREFIX/atomx_build.env"
echo "SLATE ($GPU) installed in $PREFIX"
