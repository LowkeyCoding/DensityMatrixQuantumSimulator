#include <dmqs/dmqs.hpp>
#include <dmqs/channels.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
int main(int argc, char const *argv[])
{
    const double t1[] =  {20,20};
    const double t2[] = {18,18};

    cx_mat rho = BinaryStringToDensityMatrix("0000");
    cout << rho;
    cout << "\n";
    rho = ApplyGate(rho, GH, 0);
    cout << rho;
    cout << "Happ\n";
    rho = ApplyCGate(rho, GX, 0, 1);
    cout << rho;       
    cout << "Cxapp\n";
    rho = apply_channel(rho, bit_flip_ops(0.01),4);
    cout << "Noise\n";
    cout << rho;
    cout << "\n";
    rho = ApplyCGate(rho, GX, 0, 1);
    cout << "Cxapp\n";
    cout << rho;
    rho = ApplyGate(rho, GH, 0);
    cout << "Happ\n";
    cout << rho;
    return 0;
}
