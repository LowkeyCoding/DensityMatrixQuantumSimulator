#include "dmqs.h"
#include "uppaal.h"
#include "gates.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"


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

        CHECK(DensityMatrixApproxEq(t1, c10,0));
        CHECK(DensityMatrixApproxEq(t2, c01,0));
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
}