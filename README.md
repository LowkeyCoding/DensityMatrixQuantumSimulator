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
UPPAAL definitions for a 2 qubit system. To scale array size to `N` qubit system use the following equation `size = pow(2, 2*N+1)`. Density matrix sizes are given in the number of qubits in the system.
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
    // Takes a binary string and overwrites rho with a density matrix encoding that state.
    // Supported basis states are |0⟩ and |1⟩. 
    void UInitBinState(double& rho[32], int rho_size, const string& bin);
    // Applies a gate to rho.
    void UApplyGate(double& rho[32], int rho_size, int gate, int target);
    // Applies a controlled gate to rho.
    void UApplyCGate(double& rho[32], int rho_size, int gate, int target, int control);
    // Applies a unitary U to rho given that u_size == rho_size.
    void UApplyUnitary(double& rho[32], int rho_size, double& U[32], int u_size);
    // Generates a new density matrix prho from rho containing the state of the target qubits.
    void UPartialTrace(double& rho[32], int rho_size, double& prho[8], int prho_size, int& targets[1], int targets_size);
    // Projects a basis state on to the target qubit in rho.
    void UBasisProjection(double& rho[32], int rho_size, int target, int base);
    // Projects basis state on the target qubits in rho. 
    void UBasisProjections(double& rho[32], int rho_size, int& targets[1], int targets_size, int state);
    // Returns an int containing the result of the measurement of the target qubits in the density matrix rho. 
    // If it is needed to project the result on to rho use UBasisProjections.
    int UPartialMeasure(double& rho[32], int rho_size, int& targets[1], int targets_size, double r);
    // Returns an int containing the result of the measurement of all the qubits in the density matrix rho. 
    // If it is needed to project the result on to rho use UBasisProjections with all qubits as targets.
    int UMeasureAll(double& rho[32], int rho_size, double random_value);
    // Applies amplitude dampening and dephasing channel on rho as seen in: 10.1098/rspa.2008.0439
    void UAmplitudeDampeningAndDephasing(double& rho[32], int rho_size, double& T1[2], double& T2[2], double t);
};
```