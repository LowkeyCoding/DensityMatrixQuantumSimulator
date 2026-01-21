#include <uppaal/uppaal.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using std::cout, std::endl;


int main() {
    double rho[32] = {0};
    int rho_size = 2;

    InitBinState(rho, 2, "00");
    ApplyGate(rho, rho_size, GH, 0);
    ApplyCGate(rho, rho_size, GX, 0, 1);
    ApplyChannel(rho, rho_size, AMPLITUDE_DAMPING, 0.99);

    cout << "Final density matrix:" << endl;
    for (int i = 0; i < 32; i++) {
        cout << rho[i] << " ";
    }
    cout << endl;
    return 0;
}
