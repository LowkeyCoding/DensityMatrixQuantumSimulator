#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int main() {
    const int qc = 5;
    double* test = reinterpret_cast<double*>(calloc(131072, sizeof(double)));
    double T1[qc] = {20};
    double T2[qc] = {18};
    for (int i = 0; i < 10000; i++) {
        UInitBinState(test, qc, "10010");
        UApplyGate(test, qc, GH, 0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyCGate(test, qc, GX, 0, 1);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyGate(test, qc, GX, 0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyGate(test, qc, GZ, 0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyCGate(test, qc, GX, 0, 1);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyGate(test, qc, GH, 0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UMeasureAll(test, qc, rand() / RAND_MAX);
        std:: cout << "Iteration " << i << " complete." << std::endl;
    }
}
