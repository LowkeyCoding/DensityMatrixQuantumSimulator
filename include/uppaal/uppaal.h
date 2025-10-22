#include <dmqs/dmqs.hpp>
#include <cstdlib>
#ifndef UPPAAL_H
#define UPPAAL_H

cx_mat ToMatrix(double* matrix, int size);
void FromMatrix(cx_mat matrix, double* ret, int size);

extern "C" void UInitBinState(double* rho, int rho_size, const char* state);
extern "C" int UMeasureAll(double* rho, int rho_size, double random_value);
extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target);
extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int target, int control);
extern "C" void UApplyMGate(double* rho, int rho_size, int target, double random);
extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t);

// New Format
// Assume all functions retuns the id of a given matrix operation
/*
int init_state(double* state, size_t size);
int init_state(const char* state, size_t size);
int apply_gate(int rho, int gate, int target);
int apply_noise_addph(int rho, double* t1, double* t2, double t);
int meassure_all(int rho);
int meassure(int rho, int* qubits);
int trace(int rho);
bool is_pure(int rho);
*/
#endif // UPPAAL_H