#include "dmqs.h"
#ifdef __cplusplus
extern "C" { // tells C++ compiler to use C symbol name mangling (C compiler ignores)
#endif // __cplusplus
void UMeasureAll(double* rho, int rho_size, double* random_values, int32_t* res, int32_t smaple_count);
void UApplyGate(double* rho, int rho_size, int32_t gate, int32_t target);
void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t);
cx_mat ToMatrix(double* matrix, int32_t size);
void FromMatrix(cx_mat matrix, double* ret, int size);
#ifdef __cplusplus
} // end of "C" symbol name mangling
#endif // __cplusplus