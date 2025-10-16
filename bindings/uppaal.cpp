#include "uppaal.h"
#include <iostream>
#include <fstream>
extern "C" void UInitBinState(double* rho, int rho_size, int state[]) {
    //assert(size_t(rho_size) == strlen(state));
    std::string str;
    for(int i =0;i<rho_size;++i)
      str.push_back('0'+state[i]);

    cx_mat res = BinaryStringToDensityMatrix(str);
    FromMatrix(res, rho, rho_size);
}

extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target) {
    cx_mat res = ApplyGate(ToMatrix(rho,rho_size), static_cast<u_gate>(gate), target);
    FromMatrix(res, rho, rho_size);
}

extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t){
    auto in_mat = ToMatrix(rho,rho_size);
    cx_mat res = ApplyAmplitudeDampeningAndDephasing(in_mat, T1,T2, t);
    ofstream myfile;
    myfile.open("/home/lokew/Documents/code/DensityMatrixQuantumSimulator/log.txt");
    myfile << "Input matrix: \n";
    myfile << in_mat << "\n";
    myfile << "T1:" << T1[0] << "," << T1[2] << "\n";
    myfile << "T2:" << T2[0] << "," << T2[2] << "\n";
    myfile << "t:" << t << "\n";
    myfile << "Noisy matrix:\n";
    myfile << res << "\n";
    myfile.close();
    FromMatrix(res, rho, rho_size);
}

extern "C" void UMeasureAllS(double* rho, int rho_size, double* random_values, int* res, int smaple_count) {
    vector<double> rv = vector<double>(random_values, (random_values + sizeof(double)*smaple_count));    
    for(int i = 0; i < smaple_count; i++) {
        res[i] = Sample(ToMatrix(rho,rho_size),rv[i]);
    }
}

extern "C" int UMeasureAll(double* rho, int rho_size, double random_value) {
    return Sample(ToMatrix(rho,rho_size),random_value);
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
