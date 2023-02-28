# marker-parent-set

## Build with CMake

```console
foo@bar:~/mps$ cmake -S . -B build
foo@bar:~/mps$ cmake --build build
```

## Run tests

```console
foo@bar:~/mps$ cd build && ctest
```

## Benchmark

Our reference is Plink's ld functionality. This is run with

```console
foo@bar:~/$ time plink --bfile $BEDNAME --r --matrix --out $OUTFILE
```

and compare that to

```console
foo@bar:~/ time ./mps mcorrp PREPDIR CHR
```

## Common Errors

`no kernel image is available for execution on the device`

Can be caused by the chosen device having a lower #[GPU Compute Capability](https://developer.nvidia.com/cuda-gpus#collapseOne) than the one MPS was compiled for. The Compute Capability targeted by the build is specified in the top level `CMakeLists.txt`. If there are multiple devices on the machine and only a subset of them have the appropriate Compute Capability, you can choose one by setting `CUDA_VISIBLE_DEVICES=X` where X is the index of the device. The device list can be inspected with `nvidia-smi -L`.

