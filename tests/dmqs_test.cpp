#include "../src/dmqs.hpp"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../libs/doctest.h"


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
    SUBCASE("Apply gate using string") {
        // TODO: Add more tests
        cx_mat t1 = ApplyGate(c00, GX, 0);
        cx_mat t2 = ApplyGate(c00, GX, 1);
        cx_mat t3 = ApplyGate(c00, GCX, 0);
        cx_mat t4 = ApplyGate(c10, GCX, 0);
        cx_mat t5 = ApplyGate(c11, GCX, 0);

        CHECK(DensityMatrixApproxEq(t1, c10,0));
        CHECK(DensityMatrixApproxEq(t2, c01,0));
        CHECK(DensityMatrixApproxEq(t3, c00,0));
        CHECK(DensityMatrixApproxEq(t4, c11,0));
        CHECK(DensityMatrixApproxEq(t5, c10,0));
    }
    

    CHECK(DensityMatrixApproxEq(m00,c00,0));
    CHECK(DensityMatrixApproxEq(m01,c01,0));
    CHECK(DensityMatrixApproxEq(m10,c10,0));
    CHECK(DensityMatrixApproxEq(m11,c11,0));
}

TEST_CASE("Rotation Gates") {
    cx_mat rz = RZ(80);
    cx_mat rz2 = RZ(80.0001);
    CHECK(DensityMatrixApproxEq(rz, rz, 0.0001));
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
        auto CZ12_T1 = ApplyGateToDensityMatrix(states["10"], CG(Y(), 1, 2));
        auto CZ12_F1 = ApplyGateToDensityMatrix(states["10"], CG(Z(), 1, 2));
        cout << states["10"] << endl;
        cout << CZ12_T1 << endl;
        cout << CZ12_F1 << endl;
    }
    SUBCASE("q2 > q1") {
        auto C21_T1 = ApplyGateToDensityMatrix(states["11"], CG(X(), 2, 1));
        auto C21_F1 = ApplyGateToDensityMatrix(states["10"], CG(X(), 2, 1));
        
        auto C21_T0 = ApplyGateToDensityMatrix(states["01"], CG(X(), 2, 1));
        auto C21_F0 = ApplyGateToDensityMatrix(states["00"], CG(X(), 2, 1));

        auto C31_T1 = ApplyGateToDensityMatrix(states["111"], CG(X(), 3, 1));
        auto C31_F1 = ApplyGateToDensityMatrix(states["110"], CG(X(), 3, 1));

        auto C31_T0 = ApplyGateToDensityMatrix(states["011"], CG(X(), 3, 1));
        auto C31_F0 = ApplyGateToDensityMatrix(states["010"], CG(X(), 3, 1));

        // Should create a function that given the state, gate and the first qubit to attach to will resize the gate such that both matrixes are of the same size.
        auto C32_T1 = ApplyGateToDensityMatrix(states["011"], kron(Id(),CG(X(), 3, 2)));
        auto C32_F1 = ApplyGateToDensityMatrix(states["010"], kron(Id(),CG(X(), 3, 2)));
        
        // Is CX 1,2 still CX
        CHECK(DensityMatrixApproxEq(CG(X(),1,2), CX(), 0.0));
        
        CHECK(DensityMatrixApproxEq(C21_T1, states["01"], 0.0));
        CHECK(DensityMatrixApproxEq(C21_F1, states["10"], 0.0));

        CHECK(DensityMatrixApproxEq(C21_T0, states["11"], 0.0));
        CHECK(DensityMatrixApproxEq(C21_F0, states["00"], 0.0));

        CHECK(DensityMatrixApproxEq(C31_T1, states["011"], 0.0));
        CHECK(DensityMatrixApproxEq(C31_F1, states["110"], 0.0));

        CHECK(DensityMatrixApproxEq(C31_T0, states["111"], 0.0));
        CHECK(DensityMatrixApproxEq(C31_F0, states["010"], 0.0));

        CHECK(DensityMatrixApproxEq(C32_T1, states["001"], 0.0));
        CHECK(DensityMatrixApproxEq(C32_F1, states["010"], 0.0));
    }
}
TEST_CASE("Rearrange bits"){
    vector<size_t> t12 ={1,2};
    vector<size_t> t03 ={0,3};
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
    CHECK(DensityMatrixApproxEq(PartialTrace(rho,t12),c00,0));
    CHECK(DensityMatrixApproxEq(PartialTrace(rho,t03),c11,0));
}

TEST_CASE("Sample Bell State"){
    cx_mat rho = BinaryStringToDensityMatrix("00");
    rho = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(H(),0,2));
    rho = ApplyGateToDensityMatrix(rho, CX());
    vector<double> random = {0.31139488122381065,0.706346140299032,0.0811185803790464,0.8431768096973808,0.9705136993318086,0.8775436953654865,0.6858954365356232,0.18346957566198463,0.2716272931656005,0.8928967556999657,};
    size_t counts[] = {0,0,0,0};
    size_t pre_counts[] = {4,0,0,6};
    for (size_t v : Sample(rho,random)) {
       counts[v] += 1;
    }
    for(size_t i = 0; i < 4; i++) {
        CHECK(pre_counts[i] == counts[i]);
    }
    rho = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(H(),0,2));
    SUBCASE("Pure and Mixed States") {
        INFO("Trace:", trace(rho*rho));
        double T1[] = {18,18};
        double T2[] = {20,20};
        CHECK(IsPure(ApplyAmplitudeDampeningAndDephasing(rho, T1, T2, 100), 0.000000000000001));
    }
}

