#include "dmqs.h"
#ifdef __cplusplus
extern "C" { // tells C++ compiler to use C symbol name mangling (C compiler ignores)
#endif // __cplusplus
void UInitBinState(double* rho, int rho_size, const char* state);
void UMeasureAll(double* rho, int rho_size, double* random_values, int* res, int smaple_count);
void UApplyGate(double* rho, int rho_size, int gate, int target);
void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t);
cx_mat ToMatrix(double* matrix, int size);
void FromMatrix(cx_mat matrix, double* ret, int size);
#ifdef __cplusplus
} // end of "C" symbol name mangling
#endif // __cplusplus

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