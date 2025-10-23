#include <dmqs/dmqs.hpp>
#include <cstdlib>

/*

void UInitBinState(double& rho[32], int rho_size, const string& bin);
void UApplyGate(double& rho[32], int rho_size, int gate, int target);
void UApplyCGate(double& rho[32], int rho_size, int gate, int target, int control);
void UAmplitudeDampeningAndDephasing(double& rho[32], int rho_size, double& T1[2], double& T2[2], double t);
int UMeasureAll(double& rho[32], int rho_size, double random_value);

*/

#ifndef UPPAAL_H
#define UPPAAL_H

extern "C" void UInitBinState(double* rho, int rho_size, const char* state);
extern "C" int UMeasureAll(double* rho, int rho_size, double random_value);
extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target);
extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int target, int control);
extern "C" void UApplyMGate(double* rho, int rho_size, int target, double random);
extern "C" void UApplyUnitary(double* rho, int rho_size, double* U, int u_size);
extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size, double* T1, double* T2, double t);
extern "C" void UPartialTrace(double* rho, int rho_size, double* prho, int prho_size, int* targets, int targets_size);

#endif // UPPAAL_H

