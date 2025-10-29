#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    const int qc = 14;
    double* test = (double*)calloc(536870912, sizeof(double));
    double T1[qc] = {20};
    double T2[qc] = {18};
    for(int i = 0; i < 1; i++) {
        UInitBinState(test, qc, "10001001010110");
        UApplyGate(test, qc, GH, 0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyCGate(test, qc, GX, 0, 1);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyGate(test, qc, GX,0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyGate(test, qc, GZ,0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyCGate(test, qc, GX, 0, 1);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UApplyGate(test, qc, GH, 0);
        UAmplitudeDampeningAndDephasing(test, qc, T1, T2, 1);
        UMeasureAll(test, qc, (double)rand() / (double)RAND_MAX);
    }
}