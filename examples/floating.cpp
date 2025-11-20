#include <iostream>
#include <dmqs/dmqs.hpp>

int main() {
    cx_mat result = BinaryStringToDensityMatrix("0");
    int loop_count = 1 << 10;
    std::cout << "Result: " << result << std::endl;
    int interval = 1 << 15;
    double ad = 1e-5;
    cx_mat t = BinaryStringToDensityMatrix("+");
    for (int i = 0; i < loop_count-1; i++) {
        result = ApplyGate(result, GH, 0);
        if ((i % interval) == 0) {
            std::cout << "Result: " << result << " after " << i
                      << " applications of H" << std::endl;
            std::cout << "Difference: " << result - t << std::endl;
        }
        if ((i % 2 == 0) && !mat_eq(result, t, ad)) {
            std::cout << "Error after " << i
                      << " applications of H" << std::endl;
            std::cout << "Result: " << result << std::endl;
            std::cout << "Expected: " << t << std::endl;
            std::cout << "Difference: " << result - t << std::endl;
            std::cout << "Allowed difference: " << ad << std::endl;
            return 1;
        }
    }
}
