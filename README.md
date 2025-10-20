# DensitryMatrixQuantumSimulator

## Build

Configure:
```shell
cmake -B build-release -DCMAKE_BUILD_TYPE=RelWithDebInfo
```

Compile:
```shell
cmake --build build-release
```

Test:
```shell
ctest --test-dir build-release --output-on-failure
```
