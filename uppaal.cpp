#include "uppal.h"

map<string, cx_mat> circuits;

void NewCircuit(string name, double* matrix, int size) {
    circuits.insert({name, ToMatrix(matrix, size)});
}

void UApplyGate(string circuit, string gate, int target) {
    auto elm = circuits.find(circuit);
    assert(elm != circuits.end());
    cx_mat rho = elm->second;
    rho = ApplyGate(rho, gate, target);
}

void UAmplitudeDampeningAndDephasing(string circuit, double* T1, double* T2, double t){
    auto elm = circuits.find(circuit);
    assert(elm != circuits.end());
    cx_mat rho = elm->second;
    rho = ApplyAmplitudeDampeningAndDephasing(rho, T1,T2, t);
}

void UMeasureAll(string circuit, double* random_values, int* res, int size) {
    auto elm = circuits.find(circuit);
    assert(elm != circuits.end());
    cx_mat rho = elm->second;
    vector<double> rv = vector<double>(random_values, (random_values + sizeof(double)*size));

    vector<int32_t> vres = Sample(rho, rv);
    for(size_t i = 0; i < vres.size(); i++) {
        res[i] = vres[i];
    }
}

cx_mat ToMatrix(double* matrix, int size){
    cx_double test[size*size];
    for (int i = 0; i < size*size; i++)
    {
        test[i] = cx_double(matrix[i*2], matrix[i*2+1]); 
    }
    
    cx_mat test2 = cx_mat(test,size,size, true,true);

    return test2;
}

void FromMatrix(cx_mat matrix, int size, double* ret){
    cx_double* result = matrix.memptr();
    for (int i = 0; i < size*size; i++)
    {
        ret[2*i] = result[i].real();
        ret[2*i+1] = result[i].imag(); 
    }
}
