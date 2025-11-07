#include <dmqs/dmqs.hpp>
#include <vector>
#include <string>
#include "doctest/doctest.h"

#define EXACT 0.0
#define DEC14 1e-14

TEST_CASE("Generating 1 qubit state from basis string") {
    auto q0 = BinaryStringToDensityMatrix("0");
    auto q1 = BinaryStringToDensityMatrix("1");
    auto qp = BinaryStringToDensityMatrix("+");
    auto qm = BinaryStringToDensityMatrix("-");

    cx_double d0 = cx_double(0, 0);
    cx_double d1 = cx_double(1, 0);
    cx_double dhp = cx_double(0.5, 0.0);
    cx_double dhm = cx_double(-0.5, -0.0);

    cx_mat::fixed<2, 2> c0 = {d1, d0,
                             d0, d0};
    cx_mat::fixed<2, 2> c1 = {d0, d0,
                             d0, d1};
    cx_mat::fixed<2, 2> cp = {dhp, dhp,
                             dhp, dhp};
    cx_mat::fixed<2, 2> cm = {dhp, dhm,
                             dhm, dhp};

    INFO("q0: \n", q0);
    INFO("q1: \n", q1);
    INFO("qp: \n", qp, "\ncp: \n", cp);
    INFO("qm: \n", qm, "\ncm: \n", cm);
    CHECK(mat_eq(q0, c0, EXACT));
    CHECK(mat_eq(q1, c1, EXACT));
    // As close to exact as double allows for.
    CHECK(mat_eq(qp, cp, 1e-15));
    CHECK(mat_eq(qm, cm, 1e-15));
}

TEST_CASE("Generating 2 qubit state from binary string") {
    cx_mat m00 = BinaryStringToDensityMatrix("00");
    cx_mat m01 = BinaryStringToDensityMatrix("01");
    cx_mat m10 = BinaryStringToDensityMatrix("10");
    cx_mat m11 = BinaryStringToDensityMatrix("11");
    cx_double d0 = cx_double(0, 0);
    cx_double d1 = cx_double(1, 0);

    cx_mat::fixed<4, 4> c00 = {
        d1, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0
    };
    cx_mat::fixed<4, 4> c01 = {
        d0, d0, d0, d0,
        d0, d1, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0
    };
    cx_mat::fixed<4, 4> c10 = {
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d1, d0,
        d0, d0, d0, d0
    };
    cx_mat::fixed<4, 4> c11 = {
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d1
    };
    SUBCASE("Apply gate using u_gate") {
        cx_mat t1 = ApplyGate(c00, GX, 0);
        cx_mat t2 = ApplyGate(c00, GX, 1);
        cx_mat t3 = ApplyCGate(c00, GX, 0, 1);
        cx_mat t4 = ApplyCGate(c10, GX, 0, 1);
        cx_mat t5 = ApplyCGate(c11, GX, 0, 1);

        CHECK(mat_eq(t1, c10, EXACT));
        CHECK(mat_eq(t2, c01, EXACT));
        CHECK(mat_eq(t3, c00, EXACT));
        CHECK(mat_eq(t4, c11, EXACT));
        CHECK(mat_eq(t5, c10, EXACT));
    }

    CHECK(mat_eq(m00, c00, EXACT));
    CHECK(mat_eq(m01, c01, EXACT));
    CHECK(mat_eq(m10, c10, EXACT));
    CHECK(mat_eq(m11, c11, EXACT));
}

TEST_CASE("Rotation Gates") {
    cx_mat rz = RZ(80);
    cx_mat rz2 = RZ(80);
    CHECK(mat_eq(rz, rz, EXACT));
}

TEST_CASE("Controlled Gates Test") {
    map<string, cx_mat> states = {};
    // Generate all binary strings of size 2 and 3.
    for (int n = 2; n < 4; n++) {
        for (int i = 0; i < (1 << n); i++) {
            string s;

            // build string from bits of i
            for (int j = n - 1; j >= 0; --j)
                s += ((i >> j) & 1) ? '1' : '0';
            states.insert({s, BinaryStringToDensityMatrix(s)});
        }
    }

    SUBCASE("Generated CX is CX") {
        CHECK(mat_eq(CG(X(), 0, 1), CX(), EXACT));
    }

    SUBCASE("Target > Control") {
        auto C01_T1 = ApplyGateToDensityMatrix(states["11"], CG(X(), 0, 1));
        auto C01_F1 = ApplyGateToDensityMatrix(states["01"], CG(X(), 0, 1));

        auto C01_T0 = ApplyGateToDensityMatrix(states["10"], CG(X(), 0, 1));
        auto C01_F0 = ApplyGateToDensityMatrix(states["00"], CG(X(), 0, 1));

        auto C02_T1 = ApplyGateToDensityMatrix(states["111"], CG(X(), 0, 2));
        auto C02_F1 = ApplyGateToDensityMatrix(states["011"], CG(X(), 0, 2));

        auto C02_T0 = ApplyGateToDensityMatrix(states["110"], CG(X(), 0, 2));
        auto C02_F0 = ApplyGateToDensityMatrix(states["010"], CG(X(), 0, 2));

        /*  Should create a function that given the state, gate and the 
            first qubit to attach to will resize the gate such that both 
            matrixes are of the same size. */
        auto C12_T1 = ApplyGateToDensityMatrix(
            states["011"],
            kron(Id(), CG(X(), 1, 2)));
        auto C12_F1 = ApplyGateToDensityMatrix(
            states["001"],
            kron(Id(), CG(X(), 1, 2)));

        CHECK(mat_eq(C01_T1, states["10"], EXACT));
        CHECK(mat_eq(C01_F1, states["01"], EXACT));

        CHECK(mat_eq(C01_T1, ApplyCGate(states["11"], GX, 0, 1), EXACT));
        CHECK(mat_eq(C01_F1, ApplyCGate(states["01"], GX, 0, 1), EXACT));

        CHECK(mat_eq(C01_T0, states["11"], EXACT));
        CHECK(mat_eq(C01_F0, states["00"], EXACT));

        CHECK(mat_eq(C01_T0, ApplyCGate(states["10"], GX, 0, 1), EXACT));
        CHECK(mat_eq(C01_F0, ApplyCGate(states["00"], GX, 0, 1), EXACT));

        CHECK(mat_eq(C02_T1, states["110"], EXACT));
        CHECK(mat_eq(C02_F1, states["011"], EXACT));

        CHECK(mat_eq(C02_T1, ApplyCGate(states["111"], GX, 0, 2), EXACT));
        CHECK(mat_eq(C02_F1, ApplyCGate(states["011"], GX, 0, 2), EXACT));

        CHECK(mat_eq(C02_T0, states["111"], EXACT));
        CHECK(mat_eq(C02_F0, states["010"], EXACT));

        CHECK(mat_eq(C02_T0, ApplyCGate(states["110"], GX, 0, 2), EXACT));
        CHECK(mat_eq(C02_F0, ApplyCGate(states["010"], GX, 0, 2), EXACT));

        CHECK(mat_eq(C12_T1, states["010"], EXACT));
        CHECK(mat_eq(C12_F1, states["001"], EXACT));

        CHECK(mat_eq(C12_T1, ApplyCGate(states["011"], GX, 1, 2), EXACT));
        CHECK(mat_eq(C12_F1, ApplyCGate(states["001"], GX, 1, 2), EXACT));
    }

    SUBCASE("Control > Target") {
        auto C10_T1 = ApplyGateToDensityMatrix(states["11"], CG(X(), 1, 0));
        auto C10_F1 = ApplyGateToDensityMatrix(states["10"], CG(X(), 1, 0));

        auto C10_T0 = ApplyGateToDensityMatrix(states["01"], CG(X(), 1, 0));
        auto C10_F0 = ApplyGateToDensityMatrix(states["00"], CG(X(), 1, 0));

        auto C20_T1 = ApplyGateToDensityMatrix(states["111"], CG(X(), 2, 0));
        auto C20_F1 = ApplyGateToDensityMatrix(states["110"], CG(X(), 2, 0));

        auto C20_T0 = ApplyGateToDensityMatrix(states["011"], CG(X(), 2, 0));
        auto C20_F0 = ApplyGateToDensityMatrix(states["010"], CG(X(), 2, 0));

        /*  Should create a function that given the state, gate and the first
            qubit to attach to will resize the gate such that both matrixes are 
            of the same size. */
        auto C21_T1 = ApplyGateToDensityMatrix(
            states["011"],
            kron(Id(), CG(X(), 2, 1)));

        auto C21_F1 = ApplyGateToDensityMatrix(
            states["010"],
            kron(Id(), CG(X(), 2, 1)));

        CHECK(mat_eq(C10_T1, states["01"], EXACT));
        CHECK(mat_eq(C10_F1, states["10"], EXACT));

        CHECK(mat_eq(C10_T1, ApplyCGate(states["11"], GX, 1, 0), EXACT));
        CHECK(mat_eq(C10_F1, ApplyCGate(states["10"], GX, 1, 0), EXACT));

        CHECK(mat_eq(C10_T0, states["11"], EXACT));
        CHECK(mat_eq(C10_F0, states["00"], EXACT));

        CHECK(mat_eq(C10_T0, ApplyCGate(states["01"], GX, 1, 0), EXACT));
        CHECK(mat_eq(C10_F0, ApplyCGate(states["00"], GX, 1, 0), EXACT));

        CHECK(mat_eq(C20_T1, states["011"], EXACT));
        CHECK(mat_eq(C20_F1, states["110"], EXACT));

        CHECK(mat_eq(C20_T1, ApplyCGate(states["111"], GX, 2, 0), EXACT));
        CHECK(mat_eq(C20_F1, ApplyCGate(states["110"], GX, 2, 0), EXACT));

        CHECK(mat_eq(C20_T0, states["111"], EXACT));
        CHECK(mat_eq(C20_F0, states["010"], EXACT));

        CHECK(mat_eq(C20_T0, ApplyCGate(states["011"], GX, 2, 0), EXACT));
        CHECK(mat_eq(C20_F0, ApplyCGate(states["010"], GX, 2, 0), EXACT));

        CHECK(mat_eq(C21_T1, states["001"], EXACT));
        CHECK(mat_eq(C21_F1, states["010"], EXACT));

        CHECK(mat_eq(C21_T1, ApplyCGate(states["011"], GX, 2, 1), EXACT));
        CHECK(mat_eq(C21_F1, ApplyCGate(states["010"], GX, 2, 1), EXACT));
    }
}

TEST_CASE("Rearrange bits") {
    cx_double d0 = cx_double(0, 0);
    cx_double d1 = cx_double(1, 0);
    cx_mat::fixed<4, 4> c00 = {
        d1, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0
    };

    cx_mat::fixed<4, 4> c11 = {
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d0,
        d0, d0, d0, d1
    };

    cx_mat rho = BinaryStringToDensityMatrix("1001");
    CHECK(mat_eq(PartialTrace(rho, {1, 2}), c00, EXACT));
    CHECK(mat_eq(PartialTrace(rho, {0, 3}), c11, EXACT));
    CHECK(mat_eq(PartialTrace(rho, {0}), B1(), EXACT));
    CHECK(mat_eq(PartialTrace(rho, {1}), B0(), EXACT));
    CHECK(mat_eq(PartialTrace(rho, {2}), B0(), EXACT));
    CHECK(mat_eq(PartialTrace(rho, {3}), B1(), EXACT));
    cx_mat rho1 = BinaryStringToDensityMatrix("00");
    CHECK(mat_eq(PartialTrace(rho1, {0}), B0(), DEC14));
    CHECK(mat_eq(PartialTrace(rho1, {1}), B0(), DEC14));
    rho1 = ApplyGate(rho1, GX, 0);
    CHECK(mat_eq(PartialTrace(rho1, {0}), B1(), DEC14));
    CHECK(mat_eq(PartialTrace(rho1, {1}), B0(), DEC14));
}

TEST_CASE("Sample Bell State") {
    cx_mat rho = BinaryStringToDensityMatrix("00");
    rho = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(H(), 0, 2));
    rho = ApplyGateToDensityMatrix(rho, CX());
    vector<double> random = {0.31139488122381065, 0.706346140299032,
                             0.0811185803790464, 0.8431768096973808,
                             0.9705136993318086, 0.8775436953654865,
                             0.6858954365356232, 0.18346957566198463,
                             0.2716272931656005, 0.8928967556999657};
    int counts[] = {0, 0, 0, 0};
    const int pre_counts[] = {4, 0, 0, 6};
    for (int i = 0; i < 10; i++) {
        counts[Sample(rho, random[i])] +=1;
    }
    for (int i = 0; i < 4; i++) {
        CHECK(pre_counts[i] == counts[i]);
    }
    rho = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(H(), 0, 2));
    SUBCASE("Pure and Mixed States") {
        INFO("Trace:", trace(rho*rho));
        double T1[] = {18, 18};
        double T2[] = {20, 20};
        CHECK_FALSE(IsPure(
            ApplyAmplitudeDampeningAndDephasing(rho, T1, T2, 100),
            DEC14));
    }
}

TEST_CASE("Amplitude Dampening and Dephasing") {
    SUBCASE("Single Qubit |1⟩ State with Relaxation") {
        auto state = BinaryStringToDensityMatrix("1");
        double T1[] = {10.0};  // Short relaxation time
        double T2[] = {5.0};   // Short dephasing time
        double t = 1.0;

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                state, T1, T2, t);

        // After dampening, |1⟩ should partially decay to |0⟩
        // Liklyhood of |1⟩ should decrease
        double prop = abs(noisy_state(1, 1));
        CHECK(prop < 1.0);
        CHECK(prop > 0.95);

        // Trace should be preserved
        double trace_val = trace(noisy_state).real();
        CHECK(abs(trace_val - 1.0) < DEC14);
    }

    SUBCASE("Single Qubit |0⟩ State - No Relaxation") {
        auto state = BinaryStringToDensityMatrix("0");
        double T1[] = {10.0};
        double T2[] = {5.0};
        double t = 1.0;

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                state, T1, T2, t);

        CHECK(abs(noisy_state(0, 0)) > 0.95);

        double trace_val = trace(noisy_state).real();
        CHECK(abs(trace_val - 1.0) < DEC14);
    }

    SUBCASE("Two Qubit |11⟩ State with Different Noise Parameters") {
        auto state = BinaryStringToDensityMatrix("11");
        double T1[] = {10.0, 20.0};  // Different T1 for each qubit
        double T2[] = {5.0, 15.0};   // Different T2 for each qubit
        double t = 0.5;

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                state, T1, T2, t);

        // System should still be valid density matrix
        double trace_val = trace(noisy_state).real();
        CHECK(abs(trace_val - 1.0) < DEC14);

        // Check that matrix is Hermitian (rho = rho†)
        auto conjugate_transpose = adjoint(noisy_state);
        CHECK(mat_eq(noisy_state, conjugate_transpose, DEC14));
    }

    SUBCASE("Superposition State with Dephasing") {
        auto state = BinaryStringToDensityMatrix("+");  // |+⟩ state
        double initial_coherence = abs(state(0, 1));  // Should be 0.5

        double T1[] = {1000.0}; // Very long T1 (negligible amplitude dampening)
        double T2[] = {1.0};    // Short T2 (strong dephasing)
        double t = 2.0;

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                state, T1, T2, t);

        // Coherence should be reduced due to dephasing
        double final_coherence = abs(noisy_state(0, 1));
        CHECK(final_coherence < initial_coherence);

        // Trace should still be 1
        double trace_val = trace(noisy_state).real();
        CHECK(abs(trace_val - 1.0) < DEC14);
    }

    SUBCASE("No Time Evolution - State Nearly Unchanged") {
        auto state = BinaryStringToDensityMatrix("01");
        double T1[] = {1000.0, 1000.0};  // Very long T1
        double T2[] = {1000.0, 1000.0};  // Very long T2
        double t = 0.001;                 // Very short time

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                state, T1, T2, t);

        // State should be almost identical to original
        double difference = 0.0;
        for (uword i = 0; i < state.n_elem; ++i) {
            difference += abs(state(i) - noisy_state(i));
        }
        CHECK(difference < 1e-5);  // Very small change
    }

    SUBCASE("Complete Relaxation Limit") {
        auto state = BinaryStringToDensityMatrix("1");
        double T1[] = {1.0};      // Short T1
        double T2[] = {0.5};      // Short T2
        double t = 100.0;         // Very long time

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                state, T1, T2, t);

        // After complete relaxation, should be close to
        // maximally mixed state for |1⟩.
        // In amplitude dampening channel, |1⟩ decays to |0⟩
        for (uword i = 0; i < noisy_state.n_rows; i++) {
            CHECK(noisy_state(i, i) == 0.5);
        }
    }

    SUBCASE("Bell State with Noise") {
        // Create Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
        auto state = BinaryStringToDensityMatrix("00");
        state = ApplyGateToDensityMatrix(state, GateToNQubitSystem(H(), 0, 2));
        state = ApplyGateToDensityMatrix(state, CX());

        double T1[] = {50.0, 50.0};
        double T2[] = {25.0, 25.0};
        double t = 5.0;

        auto noisy_state = ApplyAmplitudeDampeningAndDephasing(
                            state, T1, T2, t);
        // Entanglement should decrease due to noise
        CHECK(abs(state(0, 0)) > abs(noisy_state(0, 0)));
        CHECK(abs(state(3, 3)) > abs(noisy_state(3, 3)));
        CHECK(abs(state(0, 3)) > abs(noisy_state(0, 3)));
        CHECK(abs(state(3, 0)) > abs(noisy_state(3, 0)));
        // Small probability of |10⟩ and |01⟩ because of decoherence
        CHECK(abs(noisy_state(1, 1)) > 0.0045);
        CHECK(abs(noisy_state(1, 1)) > 0.0045);
        // Still a valid density matrix
        double trace_val = trace(noisy_state).real();
        CHECK(abs(trace_val - 1.0) < DEC14);
    }
}

TEST_CASE("IsPure Function") {
    SUBCASE("Pure Basis States") {
        cx_mat state_0 = BinaryStringToDensityMatrix("0");
        cx_mat state_1 = BinaryStringToDensityMatrix("1");
        cx_mat state_plus = BinaryStringToDensityMatrix("+");
        cx_mat state_minus = BinaryStringToDensityMatrix("-");

        CHECK(IsPure(state_0, DEC14));
        CHECK(IsPure(state_1, DEC14));
        CHECK(IsPure(state_plus, DEC14));
        CHECK(IsPure(state_minus, DEC14));
    }

    SUBCASE("Pure Multi-Qubit States") {
        cx_mat state_00 = BinaryStringToDensityMatrix("00");
        cx_mat state_01 = BinaryStringToDensityMatrix("01");
        cx_mat state_10 = BinaryStringToDensityMatrix("10");
        cx_mat state_11 = BinaryStringToDensityMatrix("11");

        CHECK(IsPure(state_00, DEC14));
        CHECK(IsPure(state_01, DEC14));
        CHECK(IsPure(state_10, DEC14));
        CHECK(IsPure(state_11, DEC14));

        // Create Bell state - should be pure
        cx_mat bell_state = BinaryStringToDensityMatrix("00");
        bell_state = ApplyGateToDensityMatrix(
                        bell_state, GateToNQubitSystem(H(), 0, 2));
        bell_state = ApplyGateToDensityMatrix(bell_state, CX());
        CHECK(IsPure(bell_state, DEC14));
    }

    SUBCASE("Mixed States") {
        // Create maximally mixed state for 1 qubit
        cx_mat mixed_1q = 0.5 * BinaryStringToDensityMatrix("0") +
                        0.5 * BinaryStringToDensityMatrix("1");
        CHECK_FALSE(IsPure(mixed_1q, DEC14));

        // Create partially mixed state
        cx_mat partial_mixed = 0.8 * BinaryStringToDensityMatrix("0") +
                             0.2 * BinaryStringToDensityMatrix("1");
        CHECK_FALSE(IsPure(partial_mixed, DEC14));

        // Create maximally mixed state for 2 qubits
        cx_mat mixed_2q = 0.25 * BinaryStringToDensityMatrix("00") +
                        0.25 * BinaryStringToDensityMatrix("01") +
                        0.25 * BinaryStringToDensityMatrix("10") +
                        0.25 * BinaryStringToDensityMatrix("11");
        CHECK_FALSE(IsPure(mixed_2q, DEC14));
    }

    SUBCASE("Noisy States Become Mixed") {
        auto pure_state = BinaryStringToDensityMatrix("1");
        double T1[] = {10.0};
        double T2[] = {5.0};
        double t = 1.0;

        cx_mat noisy_state = ApplyAmplitudeDampeningAndDephasing(
                                pure_state, T1, T2, t);

        // After noise application, state should become mixed
        CHECK_FALSE(IsPure(noisy_state, DEC14));
    }

    SUBCASE("Tolerance Testing") {
        auto pure_state = BinaryStringToDensityMatrix("0");

        // Should pass with reasonable tolerance
        CHECK(IsPure(pure_state, 0.1));
        CHECK(IsPure(pure_state, 0.01));
        CHECK(IsPure(pure_state, 1e-10));
        CHECK(IsPure(pure_state, 1e-14));

        // Create a slightly impure state
        cx_mat nearly_pure = 0.999 * BinaryStringToDensityMatrix("0") +
                           0.001 * BinaryStringToDensityMatrix("1");

        // Should pass with loose tolerance
        CHECK(IsPure(nearly_pure, 0.01));

        // Should fail with tight tolerance
        CHECK_FALSE(IsPure(nearly_pure, 1e-10));
    }

    SUBCASE("After Measurement") {
        auto state = BinaryStringToDensityMatrix("+"); // |+⟩ state

        // Before measurement - pure
        CHECK(IsPure(state, DEC14));

        // After basis projection - still pure
        auto projected = BasisProjection(state, 0, 0);
        CHECK(IsPure(projected, DEC14));

        // After partial trace - becomes mixed if we trace out
        // part of entangled system
        auto bell_state = BinaryStringToDensityMatrix("00");
        bell_state = ApplyGateToDensityMatrix(
                        bell_state, GateToNQubitSystem(H(), 0, 2));
        bell_state = ApplyGateToDensityMatrix(bell_state, CX());

        // Trace out first qubit
        auto reduced = PartialTrace(bell_state, {0});
        // Reduced state of entangled system is mixed
        CHECK_FALSE(IsPure(reduced, DEC14));
    }
}
TEST_CASE("Gate Application Errors") {
    SUBCASE("Invalid Gate Type") {
        auto state = BinaryStringToDensityMatrix("0");
        CHECK_THROWS_WITH(
            ApplyGate(state, static_cast<u_gate>(999), 0),
            doctest::Contains("invalid gate"));
    }
}

TEST_CASE("Partial Trace Edge Cases") {
    SUBCASE("Empty Targets") {
        auto state = BinaryStringToDensityMatrix("00");
        vector<int> empty_targets;
        CHECK_THROWS(PartialTrace(state, empty_targets));
    }

    SUBCASE("All Qubits Traced") {
        auto state = BinaryStringToDensityMatrix("01");
        vector<int> all_qubits = {0, 1};
        auto result = PartialTrace(state, all_qubits);
        CHECK(abs(trace(result) - 1.0) < DEC14);
    }
}

TEST_CASE("Basis Projection Edge Cases") {
    SUBCASE("Project Already Projected State") {
        auto state = BinaryStringToDensityMatrix("0");
        state = BasisProjection(state, 0, 0); // Project to |0⟩
        auto after_projection = BasisProjection(state, 0, 0);
        CHECK(mat_eq(state, after_projection, DEC14));
    }
}

TEST_CASE("Quantum Teleportation") {
    cx_mat rho = BinaryStringToDensityMatrix("100");
    rho = ApplyGate(rho, GH, 0);
    cx_mat qpsi = PartialTrace(rho, {0});
    INFO("Input qubit:\n", qpsi);
    rho = ApplyGate(rho, GH, 1);
    rho = ApplyCGate(rho, GX, 1, 2);
    rho = ApplyCGate(rho, GX, 0, 1);
    rho = ApplyGate(rho, GH, 0);
    int sample = PartialSample(rho, {0, 1}, 0.8); // 11
    INFO("Sampled state: ", sample);
    rho = BasisProjections(rho, {0, 1}, sample);
    rho = ApplyGate(rho, GX, 2);
    rho = ApplyGate(rho, GZ, 2);
    cx_mat qtele = PartialTrace(rho, {2});
    INFO("State after teleportation: \n", rho);
    INFO("Final qubit:\n", qtele);
    CHECK(mat_eq(qpsi, qtele, DEC14));
}
