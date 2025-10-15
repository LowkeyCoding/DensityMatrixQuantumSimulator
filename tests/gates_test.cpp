#include "../src/gates/gates.hpp"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../libs/doctest.h"

#define EXACT 0.0
#define DEC14 1e-14

TEST_CASE("Unitary") {
    SUBCASE("Id") {
        auto U = Id();
        auto UD = adjoint(U);
        CHECK(DensityMatrixApproxEq(U*UD, Id(), EXACT));
    }
    SUBCASE("X") {
        auto U = X();
        auto UD = adjoint(U);
        CHECK(DensityMatrixApproxEq(U*UD, Id(), EXACT));
    }
    SUBCASE("Y") {
        auto U = Y();
        auto UD = adjoint(U);
        CHECK(DensityMatrixApproxEq(U*UD, Id(), EXACT));
    }
    SUBCASE("Z") {
        auto U = Z();
        auto UD = adjoint(U);
        auto UI = U.i();
        CHECK(DensityMatrixApproxEq(UD, UI, EXACT));
        CHECK(DensityMatrixApproxEq(U*UD, Id(), EXACT));
    }
    SUBCASE("B0") {
        auto U = B0();
        auto UD = adjoint(U);
        CHECK(DensityMatrixApproxEq(U*UD, Id(), EXACT));
    }
    SUBCASE("B1") {
        auto U = B1();
        auto UD = adjoint(U);
        CHECK(DensityMatrixApproxEq(U*UD, Id(), EXACT));
    }
    SUBCASE("CX") {
        auto U = CX();
        auto UD = adjoint(U);
        CHECK(DensityMatrixApproxEq(U*UD, Id(2), EXACT));
    }
    SUBCASE("Id(n)") {
        for(int i = 0; i < 10; i++){
            auto U = Id(i);
            auto UD = adjoint(U);
            CHECK(DensityMatrixApproxEq(U*UD, Id(i), EXACT));
        }
    }

    std::vector<double> angles = {0, 1, 45, 90, 180, 270, 360, 720};
    SUBCASE("RX(theta)") {
        for(double ok : angles){
            auto U = RX(ok);
            auto UD = adjoint(U);
            INFO("RX(",ok,") is : ", U);
            INFO("RX(",ok,")† is : ", UD);
            INFO("RX(",ok,")*RX(",ok,") is : ", U*UD);
            CHECK(approx_equal(U*UD, Id(), "absdiff", DEC14));
        }
    }
    SUBCASE("RY(theta)") {
        for(double ok : angles){
            auto U = RY(ok);
            auto UD = adjoint(U);
            INFO("RY(",ok,") is : ", U);
            INFO("RY(",ok,")† is : ", UD);
            INFO("RY(",ok,")*RY(",ok,") is : ", U*UD);
            CHECK(approx_equal(U*UD, Id(), "absdiff", DEC14));
        }
    }
    SUBCASE("RZ(theta)") {
        for(double ok : angles){
            auto U = RZ(ok);
            auto UD = adjoint(U);
            INFO("RZ(",ok,") is : ", U);
            INFO("RZ(",ok,")† is : ", UD);
            INFO("RZ(",ok,")*RZ(",ok,") is : ", U*UD);
            CHECK(approx_equal(U*UD, Id(), "absdiff", DEC14));
        }
    }
}