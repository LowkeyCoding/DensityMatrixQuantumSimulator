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

## UPPAAL
UPPAAL definitions for a 2 qubit system. To scale array size to `N` qubit system use the following equation `size = pow(2, 2*N+1)`.
```cpp
int ID = 0;
int X = 1;
int Y = 2;
int Z = 3;
int H = 4;
int CX = 5;
int CY = 6;
int CZ = 7;
int CH = 8;

import "build-release/bindings/libdmqs_uppaal.so" {
    void UInitBinState(double& rho[32], int rho_size, const string& bin);
    void UApplyGate(double& rho[32], int rho_size, int gate, int target);
    void UApplyCGate(double& rho[32], int rho_size, int gate, int target, int control);
    void UApplyMGate(double& rho[32], int rho_size, int target, double random);
    void UApplyUnitary(double& rho[32], int rho_size, double& U[32], int u_size);
    void UPartialTrace(double& rho[32], int rho_size, double& prho[8], int prho_size, int& targets[1], int targets_size);
    int UMeasureAll(double& rho[32], int rho_size, double random_value);
    void UAmplitudeDampeningAndDephasing(double& rho[32], int rho_size, double& T1[2], double& T2[2], double t);
};
```