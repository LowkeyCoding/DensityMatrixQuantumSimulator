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
        InitBinState(test, qc, "100");
        ApplyGate(test, qc, GH, 0);
        ApplyChannel(test, qc, AMPLITUDE_DAMPING, p1);
        ApplyChannel(test, qc, PHASE_DAMPING, p2);
        ApplyCGate(test, qc, GX, 0, 1);
        ApplyChannel(test, qc, AMPLITUDE_DAMPING, p1);
        ApplyChannel(test, qc, PHASE_DAMPING, p2);
        ApplyGate(test, qc, GX, 0);
        ApplyChannel(test, qc, AMPLITUDE_DAMPING, p1);
        ApplyChannel(test, qc, PHASE_DAMPING, p2);
        ApplyGate(test, qc, GZ, 0);
        ApplyChannel(test, qc, AMPLITUDE_DAMPING, p1);
        ApplyChannel(test, qc, PHASE_DAMPING, p2);
        ApplyCGate(test, qc, GX, 0, 1);
        ApplyChannel(test, qc, AMPLITUDE_DAMPING, p1);
        ApplyChannel(test, qc, PHASE_DAMPING, p2);
        ApplyGate(test, qc, GH, 0);
        ApplyChannel(test, qc, AMPLITUDE_DAMPING, p1);
        ApplyChannel(test, qc, PHASE_DAMPING, p2);
        MeasureAll(test, qc, rand() / RAND_MAX);
    }
}
