#include "uppaal.h"

extern "C" void UInitBinState(double* rho, int rho_size, const char* state) {
    assert(size_t(rho_size) == strlen(state));
    cx_mat res = BinaryStringToDensityMatrix(state);
    FromMatrix(res, rho, rho_size);
}

extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target) {
    cx_mat res = ApplyGate(ToMatrix(rho,rho_size), static_cast<u_gate>(gate), target);
    FromMatrix(res, rho, rho_size);
}

extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t){
    cx_mat res = ApplyAmplitudeDampeningAndDephasing(ToMatrix(rho,rho_size), T1,T2, t);
    FromMatrix(res, rho, rho_size);
}

extern "C" void UMeasureAll(double* rho, int rho_size, double* random_values, int* res, int smaple_count) {
    vector<double> rv = vector<double>(random_values, (random_values + sizeof(double)*smaple_count));    
    for(int i = 0; i < smaple_count; i++) {
        res[i] = Sample(ToMatrix(rho,rho_size),rv[i]);
    }
}

void FromMatrix(cx_mat matrix, double* ret, int size){
    cx_double* result = matrix.memptr();
    int mat_row = size*2;
    for (int i = 0; i < mat_row*mat_row; i++)
    {
        ret[2*i] = result[i].real();
        ret[2*i+1] = result[i].imag(); 
    }
}

cx_mat ToMatrix(double* matrix, int size){
    int mat_row = size*2;
    auto test = vector<cx_double>(mat_row*mat_row);
    for (int i = 0; i < mat_row*mat_row; i++)
    {
        test[i] = cx_double(matrix[i*2], matrix[i*2+1]); 
    }
    return cx_mat(&test[0],mat_row,mat_row, true,true);
}
