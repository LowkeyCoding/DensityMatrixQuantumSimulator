#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int main() {
    const int qc = 5;
    double* test = reinterpret_cast<double*>(calloc(131072, sizeof(double)));
    for (int i = 0; i < 10000; i++) {
        UInitBinState(test, qc, "10010");
        UApplyGate(test, qc, GH, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, 0.002);
        UApplySChannel(test, qc, PHASE_DAMPING, 0.002);
        UApplyCGate(test, qc, GX, 0, 1);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, 0.002);
        UApplySChannel(test, qc, PHASE_DAMPING, 0.002);
        UApplyGate(test, qc, GX, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, 0.002);
        UApplySChannel(test, qc, PHASE_DAMPING, 0.002);
        UApplyGate(test, qc, GZ, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, 0.002);
        UApplySChannel(test, qc, PHASE_DAMPING, 0.002);
        UApplyCGate(test, qc, GX, 0, 1);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, 0.002);
        UApplySChannel(test, qc, PHASE_DAMPING, 0.002);
        UApplyGate(test, qc, GH, 0);
        UApplySChannel(test, qc, AMPLITUDE_DAMPING, 0.002);
        UApplySChannel(test, qc, PHASE_DAMPING, 0.002);
        UMeasureAll(test, qc, rand() / RAND_MAX);
    }
}
