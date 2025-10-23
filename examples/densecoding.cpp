#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>
int main() {
    double test [32] = {0};
    double T1[2] = {20,20};
    double T2[2] = {18,18};
    for(int i = 0; i < 1000000; i++) {
        UInitBinState(test, 2, "00");
        UApplyGate(test, 2, GH, 0);
        UAmplitudeDampeningAndDephasing(test, 2, T1, T2, 1);
        UApplyCGate(test, 2, GCX, 0, 1);
        UAmplitudeDampeningAndDephasing(test, 2, T1, T2, 1);
        UApplyGate(test, 2, GX,0);
        UAmplitudeDampeningAndDephasing(test, 2, T1, T2, 1);
        UApplyGate(test, 2, GZ,0);
        UAmplitudeDampeningAndDephasing(test, 2, T1, T2, 1);
        UApplyCGate(test, 2, GCX, 0, 1);
        UAmplitudeDampeningAndDephasing(test, 2, T1, T2, 1);
        UApplyGate(test, 2, GH, 0);
        UAmplitudeDampeningAndDephasing(test, 2, T1, T2, 1);
        UMeasureAll(test, 2, (double)rand() / (double)RAND_MAX);
    }
}