#include <vector>
#include <iostream>

#include <dmqs/dmqs.hpp>
#include <dmqs/channels.hpp>

int main() {
    cx_mat rho = BinaryStringToDensityMatrix("01");
    double p = 0.5;
    kraus_ops channel = amplitude_damping_ops(1 - p);
    cx_mat res = apply_channel(rho, channel);
    cx_mat expected = BinaryStringToDensityMatrix("00");
    std::cout << "Kraus Operators:\n";
    for (auto K : channel) {
        std::cout << K << "\n";
    }
    std::cout << "Result:\n" << res << "\nExpected:\n" << expected << std::endl;
    std::cout << "Test " << (mat_eq(res, expected, 1e-14) ? "Passed" : "Failed")
              << std::endl;
    return 0;
}
