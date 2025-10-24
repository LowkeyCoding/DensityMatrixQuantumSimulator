#include <uppaal/uppaal.h>

extern "C" void UInitBinState(double* rho, int rho_size, const char* state) {
    assert(size_t(rho_size) == strlen(state));
    int mat_row = pow(2, rho_size);
    cx_mat res = BinaryStringToDensityMatrix(state);
    memcpy(rho,res.memptr(), (mat_row*mat_row*2)*sizeof(double));
}

extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target) {
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyGate(in_mat, static_cast<u_gate>(gate), target);
    memcpy(rho,res.memptr(), (mat_row*mat_row*2)*sizeof(double));
}

extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int target, int control) {
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyCGate(in_mat, static_cast<u_gate>(gate), target, control);
    memcpy(rho,res.memptr(), (mat_row*mat_row*2)*sizeof(double));
}

extern "C" void UApplyUnitary(double* rho, int rho_size, double* U, int u_size) {
    if (rho_size != u_size) {
        throw invalid_argument("Density matrix size should be equivilant to unitary. Got " + to_string(rho_size) + " and " + to_string(u_size)); 
    }
    int mat_row = pow(2, rho_size);
    cx_mat in_rho = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat in_U = cx_mat((cx_double*)U, mat_row,mat_row, false, true);
    cx_mat res = ApplyGateToDensityMatrix(in_rho, in_U);
    memcpy(rho,res.memptr(), (mat_row*mat_row*2)*sizeof(double));
}


extern "C" void UPartialTrace(double* rho, int rho_size, double* prho, int prho_size, int* targets, int targets_size) {
    if (prho_size != targets_size) {
        throw invalid_argument("The partial density matrix should have the same size as the number of targets. Got " + to_string(rho_size) + " and " + to_string(targets_size)); 
    }
    int mat_row = pow(2, rho_size);
    vector<int> t = vector<int>(targets, targets + targets_size);
    cx_mat in_rho = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = PartialTrace(in_rho, t);
    int pmat_row = pow(2,prho_size);
    memcpy(prho,res.memptr(), (pmat_row*pmat_row*2)*sizeof(double));
}

extern "C" void UBasisProjection(double* rho, int rho_size, int target, int base) {
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, true, true);
    cx_mat res = BasisProjection(in_mat, target, base);
    memcpy(rho,res.memptr(), (mat_row*mat_row*2)*sizeof(double));
}


extern "C" void UBasisProjections(double* rho, int rho_size, int* targets, int targets_size, int state) {
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    vector<int> t = vector<int>(targets, targets + targets_size);
    cx_mat res = BasisProjections(in_mat, t, state);
    memcpy(rho,res.memptr(), (pow(2,2*rho_size +1))*sizeof(double));
}

extern "C" int UPartialMeasure(double* rho, int rho_size, int* targets, int targets_size, double random){
  int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    vector<int> t = vector<int>(targets, targets + targets_size);
    return PartialSample(in_mat, t, random);
}

extern "C" int UMeasureAll(double* rho, int rho_size, double random_value) {
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    return Sample(in_mat,random_value);
}

extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t){
    int mat_row = pow(2, rho_size);
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyAmplitudeDampeningAndDephasing(in_mat, T1,T2, t);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}
