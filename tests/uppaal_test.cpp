#include <uppaal/uppaal.h>
#include "doctest/doctest.h"
#define EXACT 0.0
#define DEC14 1e-14

TEST_CASE("State intializer") {
    int qubits = 2;
    double r00[32] = {0};
    double r01[32] = {0};
    double r10[32] = {0};
    double r11[32] = {0};

    double c00[32] = {0};
    c00[0] = 1;
    double c01[32] = {0};
    c01[10] = 1;
    double c10[32] = {0};
    c10[20] = 1;
    double c11[32] = {0};
    c11[30] = 1;

    UInitBinState(r00, qubits, "00");
    CHECK(equal(begin(r00), end(r00), begin(c00)));
    UInitBinState(r01, qubits, "01");
    CHECK(equal(begin(r01), end(r01), begin(c01)));
    UInitBinState(r10, qubits, "10");
    CHECK(equal(begin(r10), end(r10), begin(c10)));
    UInitBinState(r11, qubits, "11");
    CHECK(equal(begin(r11), end(r11), begin(c11)));
    // | 0 0 >
    //   ⮤ (index 0)
    //     ⮤ (index 1)
    SUBCASE("sample states") {
        double random[20] = {0.0045, 0.7908, 0.9903, 0.1085, 0.8035, 0.6303,
                             0.2853, 0.1125, 0.3623, 0.5738, 0.9714, 0.1138,
                             0.9487, 0.3961, 0.0949, 0.0730, 0.3724, 0.5300,
                             0.9606, 0.8233};
        for (int i = 0; i < 20; i++) {
            CHECK_EQ(UMeasureAll(r00, qubits, random[i]), 0);
            CHECK_EQ(UMeasureAll(r01, qubits, random[i]), 1);
            CHECK_EQ(UMeasureAll(r10, qubits, random[i]), 2);
            CHECK_EQ(UMeasureAll(r11, qubits, random[i]), 3);
        }
    }
}

bool cmp(const double* a, const double* b, int size,  double delta) {
    for (int i = 0; i < size; i++) {
        if (abs(a[i]) - abs(b[i]) > delta) {
            return false;
        }
    }
    return true;
}

TEST_CASE("Apply Gate") {
    double cb0[8] = {0};
    UInitBinState(cb0, 1, "0");
    double cb1[8] = {0};
    UInitBinState(cb1, 1, "1");
    double b0[8] = {0};
    UInitBinState(b0, 1, "0");
    double b1[8] = {0};
    UInitBinState(b1, 1, "1");
    const double cbh[8] = {0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0};
    int size = 8;
    SUBCASE("X Gate") {
        UApplyGate(b0, 1, GX, 0);
        // 0 X == 1
        CHECK(cmp(b0, cb1, size, EXACT));
        UApplyGate(b0, 1, GX, 0);
        // 0 X X == 0
        CHECK(cmp(b0, cb0, size, EXACT));
        UApplyGate(b1, 1, GX, 0);
        // 1 X == 0
        CHECK(cmp(b1, cb0, size, EXACT));
        UApplyGate(b1, 1, GX, 0);
        // 1 X X == 1
        CHECK(cmp(b1, cb1, size, EXACT));
    }
    SUBCASE("Y Gate") {
        const double y1_expected[8] = {-1, 0, 0, 0, 0, 0, 0, 0};
        UApplyGate(b0, 1, GY, 0);
        // 0 Y == 1
        CHECK(cmp(b0, cb1, size, EXACT));
        UApplyGate(b0, 1, GY, 0);
        // 0 Y Y == 0
        CHECK(cmp(b0, cb0, size, EXACT));
        UApplyGate(b1, 1, GY, 0);
        // 1 Y == 0
        CHECK(cmp(b1, y1_expected, size, EXACT));
        UApplyGate(b1, 1, GY, 0);
        // 1 Y Y == 1
        CHECK(cmp(b1, cb1, size, EXACT));
    }
    SUBCASE("Z Gate") {
        UApplyGate(b0, 1, GZ, 0);
        // 0 Z == 1
        CHECK(cmp(b0, cb0, size, EXACT));
        UApplyGate(b0, 1, GZ, 0);
        // 0 Z Z == 0
        CHECK(cmp(b0, cb0, size, EXACT));
        UApplyGate(b1, 1, GZ, 0);
        // 1 Z == 0
        CHECK(cmp(b1, cb1, size, EXACT));
        UApplyGate(b1, 1, GZ, 0);
        // 1 Z Z == 1
        CHECK(cmp(b1, cb1, size, EXACT));
    }
    SUBCASE("H Gate") {
        UApplyGate(b0, 1, GH, 0);
        // 0 H == bh
        INFO(b0);
        CHECK(cmp(b0, cbh, size, EXACT));
        UApplyGate(b0, 1, GH, 0);
        // 0 H H == 0
        CHECK(cmp(b0, cb0, size, DEC14));
        UApplyGate(b1, 1, GH, 0);
        // 1 H == bh
        CHECK(cmp(b1, cbh, size, EXACT));
        UApplyGate(b1, 1, GH, 0);
        // 1 H H == 1
        CHECK(cmp(b1, cb1, size, DEC14));
    }
}

TEST_CASE("Apply Controlled Gate") {
    double cb00[32] = {0};
    UInitBinState(cb00, 2, "00");
    double cb11[32] = {0};
    UInitBinState(cb11, 2, "11");
    double b00[32] = {0};
    UInitBinState(b00, 2, "00");
    double b11[32] = {0};
    UInitBinState(b11, 2, "11");
    int size = 32;
    SUBCASE("CX 0,1") {
        UApplyCGate(b00, 2, GX, 0, 1);
        CHECK(cmp(b00, cb00, size, EXACT));
        UApplyGate(b00, 2, GX, 0);
        UApplyCGate(b00, 2, GX, 0, 1);
        CHECK(cmp(b00, cb11, size, EXACT));
    }
    SUBCASE("CX 1,0") {
        UApplyCGate(b00, 2, GX, 1, 0);
        CHECK(cmp(b00, cb00, size, EXACT));
        UApplyGate(b00, 2, GX, 1);
        UApplyCGate(b00, 2, GX, 1, 0);
        CHECK(cmp(b00, cb11, size, EXACT));
    }
}

TEST_CASE("Basis Projections") {
    double rho[32] = {0};
    double crho[32] = {0};
    int targets[2] = {0, 1};
    int qc = 2;
    int size = 32;
    int mem_size = 32 * sizeof(double);

    double r00[32] = {0};
    double r01[32] = {0};
    double r10[32] = {0};
    double r11[32] = {0};

    UInitBinState(r00, qc, "00");
    UInitBinState(r01, qc, "01");
    UInitBinState(r10, qc, "10");
    UInitBinState(r11, qc, "11");
    UInitBinState(rho, qc, "00");

    SUBCASE("00") {
        UApplyGate(rho, qc, GH, 0);
        UApplyGate(rho, qc, GH, 1);
        memcpy(crho, rho, mem_size);
        UBasisProjection(rho, qc, 1, 0);
        UBasisProjection(rho, qc, 0, 0);
        UBasisProjections(crho, qc, targets, 2, 0);
        INFO("rho: \n", rho, "\ncrho: \n", crho, "\nr00: \n", r00);
        CHECK(cmp(rho, r00, size, EXACT));
        CHECK(cmp(crho, r00, size, EXACT));
        CHECK(cmp(rho, crho, size, EXACT));
    }
    SUBCASE("01") {
        UApplyGate(rho, qc, GH, 0);
        UApplyGate(rho, qc, GH, 1);
        memcpy(crho, rho, mem_size);
        UBasisProjection(rho, qc, 1, 1);
        UBasisProjection(rho, qc, 0, 0);
        UBasisProjections(crho, qc, targets, 2, 1);
        INFO("rho: \n", rho, "\ncrho: \n", crho, "\nr01: \n", r01);
        CHECK(cmp(rho, r01, size, EXACT));
        CHECK(cmp(crho, r01, size, EXACT));
        CHECK(cmp(rho, crho, size, EXACT));
    }
    SUBCASE("10") {
        UApplyGate(rho, qc, GH, 0);
        UApplyGate(rho, qc, GH, 1);
        memcpy(crho, rho, mem_size);
        UBasisProjection(rho, qc, 0, 1);
        UBasisProjection(rho, qc, 1, 0);
        UBasisProjections(crho, qc, targets, 2, 2);
        INFO("rho: \n", rho, "\ncrho: \n", crho, "\nr10: \n", r10);
        CHECK(cmp(rho, r10, size, EXACT));
        CHECK(cmp(crho, r10, size, EXACT));
        CHECK(cmp(rho, crho, size, EXACT));
    }
    SUBCASE("11") {
        UApplyGate(rho, qc, GH, 0);
        UApplyGate(rho, qc, GH, 1);
        memcpy(crho, rho, mem_size);

        UBasisProjection(rho, qc, 1, 1);
        UBasisProjection(rho, qc, 0, 1);

        UBasisProjections(crho, qc, targets, 2, 3);
        INFO("rho: \n", rho, "\ncrho: \n", crho, "\nr11: \n", r11);
        CHECK(cmp(rho, r11, size, EXACT));
        CHECK(cmp(crho, r11, size, EXACT));
        CHECK(cmp(rho, crho, size, EXACT));
    }
}

TEST_CASE("Partial Sample") {
    double rho[512] = {0};  // 4 qubits
    int qc = 4;
    UInitBinState(rho, qc, "+-10");
    SUBCASE("One Target") {
        int targets[1] = {0};
        SUBCASE("0") {
            // collapse to 1
            int sample = UPartialMeasure(rho, qc, targets, 1, 0.51);
            CHECK_EQ(sample , 1);
            // collapse to 0
            sample = UPartialMeasure(rho, qc, targets, 1, 0.49);
            CHECK_EQ(sample , 0);
        }
        SUBCASE("1") {
            targets[0] = 1;
            // collapse to 1
            int sample = UPartialMeasure(rho, qc, targets, 1, 0.51);
            CHECK_EQ(sample , 1);
            // collapse to 0
            sample = UPartialMeasure(rho, qc, targets, 1, 0.49);
            CHECK_EQ(sample , 0);
        }
        SUBCASE("2") {
            targets[0] = 2;
            // collapse to 1
            int sample = UPartialMeasure(rho, qc, targets, 1, 0.51);
            CHECK_EQ(sample , 1);
            // collapse to 1
            sample = UPartialMeasure(rho, qc, targets, 1, 0.49);
            CHECK_EQ(sample , 1);
        }
        SUBCASE("3") {
            targets[0] = 3;
            // collapse to 0
            int sample = UPartialMeasure(rho, qc, targets, 1, 0.51);
            CHECK_EQ(sample , 0);
            // collapse to 0
            sample = UPartialMeasure(rho, qc, targets, 1, 0.49);
            CHECK_EQ(sample , 0);
        }
    }
    SUBCASE("Two Targets") {
        int targets[2] = {0};
        SUBCASE("0,1") {
            targets[1] = 1;
            int sample1 = UPartialMeasure(rho, qc, targets, 2, 0.0);
            int sample2 = UPartialMeasure(rho, qc, targets, 2, 0.26);
            int sample3 = UPartialMeasure(rho, qc, targets, 2, 0.51);
            int sample4 = UPartialMeasure(rho, qc, targets, 2, 0.76);
            CHECK_EQ(sample1, 0);
            CHECK_EQ(sample2, 1);
            CHECK_EQ(sample3, 2);
            CHECK_EQ(sample4, 3);
        }
        SUBCASE("2,3") {
            targets[0] = 2;
            targets[1] = 3;
            int sample1 = UPartialMeasure(rho, qc, targets, 2, 0.0);
            int sample2 = UPartialMeasure(rho, qc, targets, 2, 0.51);
            CHECK_EQ(sample1, 2);
            CHECK_EQ(sample2, 2);
        }
        SUBCASE("0,3") {
            targets[1] = 3;
            int sample1 = UPartialMeasure(rho, qc, targets, 2, 0.0);
            int sample2 = UPartialMeasure(rho, qc, targets, 2, 0.51);
            CHECK_EQ(sample1, 0);
            CHECK_EQ(sample2, 2);
        }
        SUBCASE("1,2") {
            targets[0] = 1;
            targets[1] = 2;
            int sample1 = UPartialMeasure(rho, qc, targets, 2, 0.0);
            int sample2 = UPartialMeasure(rho, qc, targets, 2, 0.51);
            CHECK_EQ(sample1, 1);
            CHECK_EQ(sample2, 3);
        }
    }
    SUBCASE("Three Targets") {
        int targets[3] = {0};

        SUBCASE("0,1,2") {
            targets[1] = 1;
            targets[2] = 2;
            int sample1 = UPartialMeasure(rho, qc, targets, 3, 0.0);
            int sample2 = UPartialMeasure(rho, qc, targets, 3, 0.26);
            int sample3 = UPartialMeasure(rho, qc, targets, 3, 0.51);
            int sample4 = UPartialMeasure(rho, qc, targets, 3, 0.76);
            CHECK_EQ(sample1, 1);
            CHECK_EQ(sample2, 3);
            CHECK_EQ(sample3, 5);
            CHECK_EQ(sample4, 7);
        }

        SUBCASE("1,2,3") {
            targets[0] = 1;
            targets[1] = 2;
            targets[2] = 3;
            int sample1 = UPartialMeasure(rho, qc, targets, 3, 0.0);
            int sample2 = UPartialMeasure(rho, qc, targets, 3, 0.51);
            CHECK_EQ(sample1, 2);
            CHECK_EQ(sample2, 6);
        }
    }
}

TEST_CASE("Partial Trace 2 qubit") {
    double cb0[8] = {0};
    UInitBinState(cb0, 1, "0");
    double cb1[8] = {0};
    UInitBinState(cb1, 1, "1");
    double r01[32] = {0};
    double rt[8] = {0};
    int targets[1] = {1};
    UInitBinState(r01, 2, "01");
    UPartialTrace(r01, 2, rt, 1, targets, 1);
    CHECK(cmp(rt, cb1, 8, EXACT));
    targets[0] = 0;
    UPartialTrace(r01, 2, rt, 1, targets, 1);
    CHECK(cmp(rt, cb0, 8, EXACT));
}

TEST_CASE("Partial Trace 3 qubit") {
    double cm1[32] = {0};
    UInitBinState(cm1, 2, "01");
    double cm2[32] = {0};
    UInitBinState(cm2, 2, "10");
    double r010[128] = {0};
    UInitBinState(r010, 3, "010");
    double rt[32] = {0};
    int targets[2] = {0, 1};
    UPartialTrace(r010, 3, rt, 2, targets, 2);
    CHECK(cmp(rt, cm1, 8, EXACT));
    targets[0] = 1;
    targets[1] = 2;
    UPartialTrace(r010, 3, rt, 2, targets, 2);
    CHECK(cmp(rt, cm2, 8, EXACT));
}

TEST_CASE("Amplitude Dampening and Dephasing") {
    SUBCASE("Single Qubit - Pure State") {
        double rho[8] = {0};
        UInitBinState(rho, 1, "1"); // |1⟩ state

        double T1[] = {10.0}; // Short relaxation time
        double T2[] = {5.0};  // Short dephasing time
        double t = 1.0;       // Time duration

        UAmplitudeDampeningAndDephasing(rho, 1, T1, T2, t);

        // After dampening, |1⟩ should partially decay to |0⟩
        // Trace should still be 1
        double trace = rho[0] + rho[6]; // Real parts of diagonal
        CHECK(abs(trace - 1.0) < DEC14);
    }

    SUBCASE("Two Qubits - Different T1/T2") {
        double rho[32] = {0};
        UInitBinState(rho, 2, "11"); // |11⟩ state

        double T1[] = {10.0, 20.0}; // Different T1 for each qubit
        double T2[] = {5.0, 15.0};  // Different T2 for each qubit
        double t = 0.5;

        UAmplitudeDampeningAndDephasing(rho, 2, T1, T2, t);

        // System should still be valid density matrix
        // Trace should be preserved
        double trace = rho[0] + rho[10] + rho[20] + rho[30];
        CHECK(abs(trace - 1.0) < DEC14);
    }

    SUBCASE("No Time Evolution") {
        double rho[32] = {0};
        double rho_original[32] = {0};
        UInitBinState(rho, 2, "01");
        UInitBinState(rho_original, 2, "01");

        double T1[] = {1000.0, 1000.0}; // Very long T1
        double T2[] = {1000.0, 1000.0}; // Very long T2
        double t = 0.001; // Very short time

        UAmplitudeDampeningAndDephasing(rho, 2, T1, T2, t);

        // State should be almost unchanged
        bool unchanged = true;
        for (int i = 0; i < 32; i++) {
            if (abs(rho[i] - rho_original[i]) > 0.01) {
                unchanged = false;
                break;
            }
        }
        CHECK(unchanged);
    }
}

TEST_CASE("Memory Management and Buffer Sizes") {
    SUBCASE("Buffer Size Validation") {
        // Test that functions handle buffer sizes correctly
        double small_buffer[4] = {0}; // Too small for 2 qubits
        CHECK_THROWS(UInitBinState(small_buffer, 2, "00"));
    }

    SUBCASE("Memory Overlap Protection") {
        double rho[32] = {0};
        UInitBinState(rho, 2, "00");

        // Apply gate should work without memory corruption
        UApplyGate(rho, 2, GX, 0);
        UApplyGate(rho, 2, GX, 1);

        // Final state should be valid
        double trace = rho[0] + rho[10] + rho[20] + rho[30];
        CHECK(abs(trace - 1.0) < 1e-14);
    }
}

TEST_CASE("Quantum Teleportation") {
    int qc = 3;
    double rho[128] = {0};
    UInitBinState(rho, qc, "100");
    UApplyGate(rho, qc, GH, 0);
    double qpsi[8] = {0};
    int target[1] = {0};
    UPartialTrace(rho, qc, qpsi, 1, target, 1);
    INFO("Input qubit:\n", qpsi);
    UApplyGate(rho, qc, GH, 1);
    UApplyCGate(rho, qc, GX, 1, 2);
    UApplyCGate(rho, qc, GX, 0, 1);
    UApplyGate(rho, qc, GH, 0);
    int measure[2] = {0, 1};
    // Sample 11
    int sample = UPartialMeasure(rho, qc, measure, 2, 0.8);
    INFO("Sampled state: ", sample);
    UBasisProjections(rho, qc, measure, 2, sample);
    UApplyGate(rho, qc, GX, 2);
    UApplyGate(rho, qc, GZ, 2);
    double qtele[8] = {0};
    target[0] = 2;
    UPartialTrace(rho, qc, qtele, 1, target, 1);
    INFO("State after teleportation: \n", rho);
    INFO("Final qubit:\n", qtele);
    CHECK(cmp(qpsi, qtele, 8, DEC14));
}
