#include <dmqs/dmqs.hpp>
#include "doctest/doctest.h"

#define EXACT 0.0
#define DEC14 1e-14

TEST_CASE("Test generating state from binary string") {
    cx_mat m00 = BinaryStringToDensityMatrix("00");
    cx_mat m01 = BinaryStringToDensityMatrix("01");
    cx_mat m10 = BinaryStringToDensityMatrix("10");
    cx_mat m11 = BinaryStringToDensityMatrix("11");
    cx_double d0 = cx_double(0,0);
    cx_double d1 = cx_double(1,0);
    
    cx_mat::fixed<4,4> c00 = {
        d1,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0
    };
    cx_mat::fixed<4,4> c01 = {
        d0,d0,d0,d0,
        d0,d1,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0
    };
    cx_mat::fixed<4,4> c10 = {
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d1,d0,
        d0,d0,d0,d0
    };
    cx_mat::fixed<4,4> c11 = {
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d1
    };
    SUBCASE("Apply gate using u_gate") {
        // TODO: Add more tests
        cx_mat t1 = ApplyGate(c00, GX, 0);
        cx_mat t2 = ApplyGate(c00, GX, 1);
        cx_mat t3 = ApplyCGate(c00, GCX, 0,1);
        cx_mat t4 = ApplyCGate(c10, GCX, 0,1);
        cx_mat t5 = ApplyCGate(c11, GCX, 0,1);

        CHECK(mat_eq(t1, c10, EXACT));
        CHECK(mat_eq(t2, c01, EXACT));
        CHECK(mat_eq(t3, c00, EXACT));
        CHECK(mat_eq(t4, c11, EXACT));
        CHECK(mat_eq(t5, c10, EXACT));
    }
    

    CHECK(mat_eq(m00,c00,EXACT));
    CHECK(mat_eq(m01,c01,EXACT));
    CHECK(mat_eq(m10,c10,EXACT));
    CHECK(mat_eq(m11,c11,EXACT));
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
    
    SUBCASE("Controled Z Gate ") {
        auto CZ12_T1 = ApplyGateToDensityMatrix(states["10"], CG(Y(), 0, 1));
        auto CZ12_F1 = ApplyGateToDensityMatrix(states["10"], CG(Z(), 0, 1));
    }
    SUBCASE("q2 > q1") {
        auto C21_T1 = ApplyGateToDensityMatrix(states["11"], CG(X(), 1, 0));
        auto C21_F1 = ApplyGateToDensityMatrix(states["10"], CG(X(), 1, 0));
        
        auto C21_T0 = ApplyGateToDensityMatrix(states["01"], CG(X(), 1, 0));
        auto C21_F0 = ApplyGateToDensityMatrix(states["00"], CG(X(), 1, 0));

        auto C31_T1 = ApplyGateToDensityMatrix(states["111"], CG(X(), 2, 0));
        auto C31_F1 = ApplyGateToDensityMatrix(states["110"], CG(X(), 2, 0));

        auto C31_T0 = ApplyGateToDensityMatrix(states["011"], CG(X(), 2, 0));
        auto C31_F0 = ApplyGateToDensityMatrix(states["010"], CG(X(), 2, 0));

        // Should create a function that given the state, gate and the first qubit to attach to will resize the gate such that both matrixes are of the same size.
        auto C32_T1 = ApplyGateToDensityMatrix(states["011"], kron(Id(),CG(X(), 2, 1)));
        auto C32_F1 = ApplyGateToDensityMatrix(states["010"], kron(Id(),CG(X(), 2, 1)));
        
        // Is CX 1,2 still CX
        CHECK(mat_eq(CG(X(),0,1), CX(), EXACT));
        
        CHECK(mat_eq(C21_T1, states["01"], EXACT));
        CHECK(mat_eq(C21_F1, states["10"], EXACT));

        CHECK(mat_eq(C21_T0, states["11"], EXACT));
        CHECK(mat_eq(C21_F0, states["00"], EXACT));

        CHECK(mat_eq(C31_T1, states["011"], EXACT));
        CHECK(mat_eq(C31_F1, states["110"], EXACT));

        CHECK(mat_eq(C31_T0, states["111"], EXACT));
        CHECK(mat_eq(C31_F0, states["010"], EXACT));

        CHECK(mat_eq(C32_T1, states["001"], EXACT));
        CHECK(mat_eq(C32_F1, states["010"], EXACT));
    }
}
TEST_CASE("Rearrange bits"){
    cx_double d0 = cx_double(0,0);
    cx_double d1 = cx_double(1,0);
    cx_mat::fixed<4,4> c00 = {
        d1,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0
    };
        
    cx_mat::fixed<4,4> c11 = {
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d0,
        d0,d0,d0,d1
    };
    cx_mat rho = BinaryStringToDensityMatrix("1001");
    CHECK(mat_eq(PartialTrace(rho,{1,2}),c00, EXACT));
    CHECK(mat_eq(PartialTrace(rho,{0,3}),c11, EXACT));
    CHECK(mat_eq(PartialTrace(rho,{0}),B1(), EXACT));
    CHECK(mat_eq(PartialTrace(rho,{1}),B0(), EXACT));
    CHECK(mat_eq(PartialTrace(rho,{2}),B0(), EXACT));
    CHECK(mat_eq(PartialTrace(rho,{3}),B1(), EXACT));
    cx_mat rho1 = BinaryStringToDensityMatrix("00");
    CHECK(mat_eq(PartialTrace(rho1, {0}), B0(), DEC14));
    CHECK(mat_eq(PartialTrace(rho1, {1}), B0(), DEC14));
    rho1 = ApplyGate(rho1, GX,0);
    CHECK(mat_eq(PartialTrace(rho1, {0}), B1(), DEC14));
    CHECK(mat_eq(PartialTrace(rho1, {1}), B0(), DEC14));
}

TEST_CASE("Sample Bell State"){
    cx_mat rho = BinaryStringToDensityMatrix("00");
    rho = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(H(),0,2));
    rho = ApplyGateToDensityMatrix(rho, CX());
    vector<double> random = {0.31139488122381065,0.706346140299032,0.0811185803790464,0.8431768096973808,0.9705136993318086,0.8775436953654865,0.6858954365356232,0.18346957566198463,0.2716272931656005,0.8928967556999657,};
    int counts[] = {0,0,0,0};
    int pre_counts[] = {4,0,0,6};
    for(int i = 0; i < 10; i++){
        counts[Sample(rho,random[i])] +=1;
    }
    for(int i = 0; i < 4; i++) {
        CHECK(pre_counts[i] == counts[i]);
    }
    rho = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(H(),0,2));
    SUBCASE("Pure and Mixed States") {
        INFO("Trace:", trace(rho*rho));
        double T1[] = {18,18};
        double T2[] = {20,20};
        CHECK(IsPure(ApplyAmplitudeDampeningAndDephasing(rho, T1, T2, 100), DEC14));
    }
}

TEST_CASE("Quantum Teleportation") {
    cx_mat rho = BinaryStringToDensityMatrix("100");
    rho = ApplyGate(rho, GH, 0);
    cx_mat qpsi = PartialTrace(rho,{0});
    INFO("Input qubit:\n", qpsi);
    rho = ApplyGate(rho, GH, 1);
    rho = ApplyCGate(rho, GCX,1,2);
    rho = ApplyCGate(rho, GCX,0,1);
    rho = ApplyGate(rho, GH,0);
    int sample = PartialSample(rho, {0,1}, 0.8); // 11
    INFO("Sampled state: ", sample);
    rho = BasisProjections(rho, {0,1}, sample);
    rho = ApplyGate(rho, GX, 2);
    rho = ApplyGate(rho, GZ, 2);
    cx_mat qtele = PartialTrace(rho,{2});
    INFO("State after teleportation: \n", rho);
    INFO("Final qubit:\n", qtele);
    CHECK(mat_eq(qpsi,qtele, DEC14));
}
