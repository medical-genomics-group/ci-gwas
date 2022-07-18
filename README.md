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
