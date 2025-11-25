#include <cstdlib>
#include <dmqs/dmqs.hpp>
#include <dmqs/channels.hpp>
#ifndef INCLUDE_UPPAAL_UPPAAL_H_
#define INCLUDE_UPPAAL_UPPAAL_H_

extern "C" void UInitBinState(double* rho, int rho_size, const char* state);
extern "C" int UMeasureAll(double* rho, int rho_size, double r);
extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target);
extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int control,
                            int target);
extern "C" void UApplyUnitary(double* rho, int rho_size, double* U, int u_size);
extern "C" void UPartialTrace(double* rho, int rho_size, double* prho,
                              int prho_size, int* targets, int targets_size);
extern "C" void UBasisProjection(double* rho, int rho_size, int target,
                                 int base);
extern "C" void UBasisProjections(double* rho, int rho_size, int* targets,
                                  int targets_size, int state);
extern "C" int UPartialMeasure(double* rho, int rho_size, int* targets,
                               int targets_size, double r);
extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size,
                                                const double* T1,
                                                const double* T2, double t);
extern "C" void UApplyChannel(double* rho, int rho_size, int channel,
                              double* probs, int probs_size);
extern "C" void UApplySChannel(double* rho, int rho_size, int channel,
                              double probs);
#endif // INCLUDE_UPPAAL_UPPAAL_H_
