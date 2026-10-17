# SLATE + BLAS++ + LAPACK++ installed by tools/build_slate.sh.
load("@rules_cc//cc:cc_library.bzl", "cc_library")

cc_library(
    name = "slate",
    srcs = glob([
        "lib/libslate.so*",
        "lib/liblapackpp.so*",
        "lib/libblaspp.so*",
    ]),
    hdrs = glob(["include/**"]),
    includes = ["include"],
    # MPI, OpenBLAS and (when built with CUDA) the CUDA runtime libraries.
    linkopts = [
        "-L/usr/lib/x86_64-linux-gnu/openmpi/lib",
        "-lmpi",
        "-l:libopenblas.so.0",
        "-L/usr/local/cuda-12.6/lib64",
        "-Wl,-rpath,/usr/local/cuda-12.6/lib64",
        "-Wl,--as-needed",
        "-lcudart",
        "-lcublas",
        "-lcusolver",
    ],
    visibility = ["//visibility:public"],
)
