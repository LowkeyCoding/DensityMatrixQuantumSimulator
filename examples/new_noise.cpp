#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using std::cout, std::endl;


int main() {
    double rho[32] = {0};
    int rho_size = 2;

    UInitBinState(rho, 2, "00");
    UApplyGate(rho, rho_size, GH, 0);
    UApplyCGate(rho, rho_size, GX, 0, 1);
    UApplyChannel(rho, rho_size, AMPLITUDE_DAMPING, new double[1]{0.99}, 1);

    cout << "Final density matrix:" << endl;
    for (int i = 0; i < 32; i++) {
        cout << rho[i] << " ";
    }
    cout << endl;
    return 0;
}
