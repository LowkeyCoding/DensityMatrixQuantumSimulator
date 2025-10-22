#include <uppaal/uppaal.h>

extern "C" void UInitBinState(double* rho, int rho_size, const char* state) {
    assert(size_t(rho_size) == strlen(state));
    cx_mat res = BinaryStringToDensityMatrix(state);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}

extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target) {
    int mat_row = rho_size*2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyGate(in_mat, static_cast<u_gate>(gate), target);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}

extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int target, int control) {
    int mat_row = rho_size*2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyCGate(in_mat, static_cast<u_gate>(gate), target, control);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}

extern "C" void UApplyMGate(double* rho, int rho_size, int target, double random_value) {
    int mat_row = rho_size*2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = MeasurementGate(in_mat, target, random_value);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}

extern "C" void UApplyUnitary(double* rho, int rho_size, double* U, int u_size) {
    if (rho_size != u_size) {
        throw invalid_argument("Density matrix size should be equivilant to unitary not " + to_string(rho_size) + " and " + to_string(u_size)); 
    }
    int mat_row = rho_size*2;
    cx_mat in_rho = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat in_U = cx_mat((cx_double*)U, mat_row,mat_row, false, true);
    cx_mat res = ApplyGateToDensityMatrix(in_rho, in_U);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}

extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t){
    int mat_row = rho_size*2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    cx_mat res = ApplyAmplitudeDampeningAndDephasing(in_mat, T1,T2, t);
    memcpy(rho,res.memptr(), (pow(2,rho_size*rho_size)*2)*sizeof(double));
}

extern "C" int UMeasureAll(double* rho, int rho_size, double random_value) {
    int mat_row = rho_size*2;
    cx_mat in_mat = cx_mat((cx_double*)rho, mat_row,mat_row, false, true);
    return Sample(in_mat,random_value);
}