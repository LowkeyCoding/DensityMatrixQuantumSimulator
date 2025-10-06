#include "dmqs.h"
void NewCircuit(string name, double* matrix, int size);
void UMeasureAll(string circuit, double* random_values, int* res, int size);
void UApplyGate(string circuit, string gate, int target);
void UAmplitudeDampeningAndDephasing(string circuit, double* T1, double* T2, double t);
cx_mat ToMatrix(double* matrix, int size);
void FromMatrix(cx_mat matrix, int size, double* ret);