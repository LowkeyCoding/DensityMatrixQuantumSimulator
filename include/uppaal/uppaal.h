#include <cstdlib>
#include <dmqs/dmqs.hpp>
#include <dmqs/channels.hpp>
#ifndef INCLUDE_UPPAAL_UPPAAL_H_
#define INCLUDE_UPPAAL_UPPAAL_H_

extern "C" void InitBinState(double* rho, int rho_size, const char* state);
extern "C" int MeasureAll(double* rho, int rho_size, double r);
extern "C" void ApplyGate(double* rho, int rho_size, int gate, int target);
extern "C" void ApplyCGate(double* rho, int rho_size, int gate, int control,
                            int target);
extern "C" void ApplyUnitary(double* rho, int rho_size, double* U, int u_size);
extern "C" void PartialTrace(double* rho, int rho_size, double* prho,
                              int prho_size, int* targets, int targets_size);
extern "C" void BasisProjection(double* rho, int rho_size, int target,
                                 int base);
extern "C" void BasisProjections(double* rho, int rho_size, int* targets,
                                  int targets_size, int state);
extern "C" int PartialMeasure(double* rho, int rho_size, int* targets,
                               int targets_size, double r);
extern "C" void AmplitudeDampeningAndDephasing(double* rho, int rho_size,
                                                const double* T1,
                                                const double* T2, double t);
extern "C" void ApplyChannel(double* rho, int rho_size, int channel,
                              double probs);
extern "C" void ApplyGAD(double* rho, int rho_size, double p, double g);
extern "C" void ResetQubit(double* rho, int rho_size, int qubit);
#endif // INCLUDE_UPPAAL_UPPAAL_H_
