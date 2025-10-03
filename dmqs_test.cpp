#include "dmqs.h"
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

    CHECK(DensityMatrixApproxEq(m00,c00,0));
    CHECK(DensityMatrixApproxEq(m01,c01,0));
    CHECK(DensityMatrixApproxEq(m10,c10,0));
    CHECK(DensityMatrixApproxEq(m11,c11,0));
}

TEST_CASE("Rotation Gates") {
    cx_mat rvg = RV(1,2,3);
    cx_mat rvg2 = RV(1,2,3.0001);
    CHECK(DensityMatrixApproxEq(rvg, rvg2, 0.0001));
}