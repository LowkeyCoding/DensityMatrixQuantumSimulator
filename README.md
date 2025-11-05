# DensitryMatrixQuantumSimulator

A high-performance C++ library for quantum system simulation using density matrix formalism, featuring noise modeling and UPPAAL integration for formal verification of quantum protocols.

## Dependencies

Before building, install the following dependencies for your systen.

### Ubuntu/Debian
```
sudo apt install build-essential cmake libopenblas-dev
```

### Arch
```
sudo pacman -S base-devel cmake openblas
```

### MacOS
```
brew install cmake openblas
```

## Build

### Configure
```shell
cmake -B build-release -DCMAKE_BUILD_TYPE=RelWithDebInfo
```

### Compile
```shell
cmake --build build-release
```

### Test
```shell
ctest --test-dir build-release --output-on-failure
```

## UPPAAL
UPPAAL definitions for a 2 qubit system. To scale array size to `N` qubit system use the following equation `size = 1 << 2*N+1`. Density matrix sizes are given in the number of qubits in the system.
```cpp
const int ID = 0;  // Identity
const int X  = 1;  // Pauli X
const int Y  = 2;  // Pauli Y
const int Z  = 3;  // Pauli Z
const int H  = 4;  // Hadamard

import "build-release/bindings/libdmqs_uppaal.so" {
    // Initialize density matrix with binary state string (e.g., "01" for |01⟩ or "+-" for |+-⟩)
    // rho_size = 1 << 2*N+1, bin = binary state string of length N
    void UInitBinState(double& rho[size], int rho_size, const string& bin);
    
    // Apply single-qubit gate to target qubit
    void UApplyGate(double& rho[size], int rho_size, int gate, int target);
    
    // Apply controlled gate (control → target)
    void UApplyCGate(double& rho[size], int rho_size, int gate, int control, int target);
    
    // Apply custom unitary matrix U (u_size must equal rho_size)
    void UApplyUnitary(double& rho[size], int rho_size, double& U[size], int u_size);
    
    // Partial trace: extract subsystem state into prho
    // prho_size = 1 << N*2-1 where N = len(targets)
    void UPartialTrace(double& rho[size], int rho_size, double& prho[partial_size], 
                      int prho_size, int& targets[target_count], int targets_size);
    
    // Project target qubit onto basis state (0 or 1)
    void UBasisProjection(double& rho[size], int rho_size, int target, int base);
    
    // Project multiple target qubits onto basis state (bitmask)
    void UBasisProjections(double& rho[size], int rho_size, int& targets[target_count], 
                          int targets_size, int state);
    
    // Measure target qubits and return result (does not collapse state)
    // Use UBasisProjections to collapse after measurement
    int UPartialMeasure(double& rho[size], int rho_size, int& targets[target_count], 
                       int targets_size, double r);
    
    // Measure all qubits and return result (does not collapse state)
    int UMeasureAll(double& rho[size], int rho_size, double random_value);
    
    // Apply amplitude damping and dephasing noise channel on rho as seen in: 10.1098/rspa.2008.0439
    void UAmplitudeDampeningAndDephasing(double& rho[32], int rho_size, double& T1[2], double& T2[2], double t);
};
```

## Develop

When the required tools are downloaded run pre-commit installation command to ensure that pre-commit hooks run before commiting. 
```
pre-commit install
```

#### Arch based
```
pacman -S pre-commit
```

#### Via Python
```
pip install pre-commit
```