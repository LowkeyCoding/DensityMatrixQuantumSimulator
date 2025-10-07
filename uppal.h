#include "dmqs.h"

void NewCircuit(int id, double* matrix, int size);
void UMeasureAll(int id, double* random_values, int* res, int size);
void UApplyGate(int id, const char* gate, int target);
void UAmplitudeDampeningAndDephasing(int id, double* T1, double* T2, double t);
cx_mat ToMatrix(double* matrix, int size);
void FromMatrix(cx_mat matrix, int size, double* ret);