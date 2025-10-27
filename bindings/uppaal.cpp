#include <uppaal/uppaal.h>

// Takes a binary string and overwrites rho with a density matrix encoding that state.
// Supported basis states are |0⟩ and |1⟩. 
extern "C" void UInitBinState(double* rho, int rho_size, const char* state) {
    assert(size_t(rho_size) == strlen(state));
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cx_mat res = BinaryStringToDensityMatrix(state);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}

// Applies a gate to rho.
extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target) {
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyGate(in_mat, static_cast<u_gate>(gate), target);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}

// Applies a controlled gate to rho.
extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int control, int target) {
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cout << "UapplyCGate: " << mat_size << endl;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyCGate(in_mat, static_cast<u_gate>(gate), control, target);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}

// Applies a unitary U to rho given that u_size == rho_size.
extern "C" void UApplyUnitary(double* rho, int rho_size, double* U, int u_size) {
    if (rho_size != u_size) {
        throw invalid_argument("Density matrix size should be equivilant to unitary. Got " + to_string(rho_size) + " and " + to_string(u_size)); 
    }
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cx_mat in_rho = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat in_U = cx_mat((cx_double*)U, mat_row,mat_row, false, true);
    cx_mat res = ApplyGateToDensityMatrix(in_rho, in_U);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}

// Generates a new density matrix prho from rho containing the state of the target qubits.
extern "C" void UPartialTrace(double* rho, int rho_size, double* prho, int prho_size, int* targets, int targets_size) {
    if (prho_size != targets_size) {
        throw invalid_argument("The partial density matrix should have the same size as the number of targets. Got " + to_string(rho_size) + " and " + to_string(targets_size)); 
    }
    int mat_row = pow(2, rho_size);
    vector<int> t = vector<int>(targets, targets + targets_size);
    cx_mat in_rho = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = PartialTrace(in_rho, t);
    int mat_size = pow(2,2*prho_size+1);
    memcpy(prho,res.memptr(), mat_size*sizeof(double));
}

// Projects a basis state on to the target qubit in rho.
extern "C" void UBasisProjection(double* rho, int rho_size, int target, int state) {
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, true, true);
    cx_mat res = BasisProjection(in_mat, target, state);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}

// Projects basis state on the target qubits in rho. 
extern "C" void UBasisProjections(double* rho, int rho_size, int* targets, int targets_size, int state) {
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    vector<int> t = vector<int>(targets, targets + targets_size);
    cx_mat res = BasisProjections(in_mat, t, state);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}

// Returns an int containing the result of the measurement of the target qubits in the density matrix rho. 
// If it is needed to project the result on to rho use UBasisProjections.
extern "C" int UPartialMeasure(double* rho, int rho_size, int* targets, int targets_size, double r){
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    vector<int> t = vector<int>(targets, targets + targets_size);
    return PartialSample(in_mat, t, r);
}

// Returns an int containing the result of the measurement of all the qubits in the density matrix rho. 
// If it is needed to project the result on to rho use UBasisProjections with all qubits as targets.
extern "C" int UMeasureAll(double* rho, int rho_size, double r) {
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    return Sample(in_mat, r);
}

// Applies amplitude dampening and dephasing channel on rho as seen in: 10.1098/rspa.2008.0439
extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t){
    int mat_row = pow(2, rho_size);
    int mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyAmplitudeDampeningAndDephasing(in_mat, T1,T2, t);
    memcpy(rho, res.memptr(), mat_size*sizeof(double));
}
