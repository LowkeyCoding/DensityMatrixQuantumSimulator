#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int main() {
    const int qc = 3; // number of qubits
    double T1 = 20;
    double T2 = 18;
    double t  = 5;
    double p1 = exp(-t/T1);
    double p2 = exp(-t/T2);
    size_t size = 1 << (2*qc+1);
    double* test = reinterpret_cast<double*>(calloc(size, sizeof(double)));
    for (int i = 0; i < 1000000; i++) {
        UInitBinState(test, qc, "100");
        UApplyGate(test, qc, GH, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, p1);
        UApplySChannel(test, qc, PHASE_DAMPING, p2);
        UApplyCGate(test, qc, GX, 0, 1);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, p1);
        UApplySChannel(test, qc, PHASE_DAMPING, p2);
        UApplyGate(test, qc, GX, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, p1);
        UApplySChannel(test, qc, PHASE_DAMPING, p2);
        UApplyGate(test, qc, GZ, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, p1);
        UApplySChannel(test, qc, PHASE_DAMPING, p2);
        UApplyCGate(test, qc, GX, 0, 1);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, p1);
        UApplySChannel(test, qc, PHASE_DAMPING, p2);
        UApplyGate(test, qc, GH, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, p1);
        UApplySChannel(test, qc, PHASE_DAMPING, p2);
        UMeasureAll(test, qc, rand() / RAND_MAX);
    }
}
