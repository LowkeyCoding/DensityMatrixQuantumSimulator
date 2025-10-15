#include "../src/gates/gates.hpp"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../libs/doctest.h"


TEST_CASE("Unitary") {
    SUBCASE("Id") {
        CHECK(DensityMatrixApproxEq(Id()*Id(), Id(), 0.0));
    }
    SUBCASE("X") {
        CHECK(DensityMatrixApproxEq(X()*X(), Id(), 0.0));
    }
    SUBCASE("Y") {
        CHECK(DensityMatrixApproxEq(Y()*Y(), Id(), 0.0));
    }
    SUBCASE("Z") {
        CHECK(DensityMatrixApproxEq(Z()*Z(), Id(), 0.0));
    }
}