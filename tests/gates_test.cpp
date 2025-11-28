#include <dmqs/gates.hpp>
#include <vector>
#include "doctest/doctest.h"


#define EXACT 0.0
#define DEC14 1e-14

TEST_CASE("Exceptions") {
    SUBCASE("Id") {
        REQUIRE_THROWS_WITH(
            Id(0),
            doctest::Contains("N must be at least 1"));
        REQUIRE_THROWS_WITH(
            Id(-1),
            doctest::Contains("N must be at least 1"));
        REQUIRE_THROWS_WITH(
            Id(-1321),
            doctest::Contains("N must be at least 1"));
    }
    SUBCASE("CG") {
        for (int i = 0; i < 10; i++) {
            REQUIRE_THROWS_WITH(
                CG(X(), i, i),
                doctest::Contains("qubit must be different"));
        }
    }
    SUBCASE("SWAP") {
        for (int i = 0; i < 10; i++) {
            REQUIRE_THROWS_WITH(
                SWAP(i, i),
                doctest::Contains("The qubits to swap"));
        }
    }
}

TEST_CASE("Unitary") {
    SUBCASE("Id") {
        auto U = Id();
        auto UD = adjoint(U);
        CHECK(mat_eq(U * UD, Id(), EXACT));
    }
    SUBCASE("X") {
        auto U = X();
        auto UD = adjoint(U);
        CHECK(mat_eq(U * UD, Id(), EXACT));
    }
    SUBCASE("Y") {
        auto U = Y();
        auto UD = adjoint(U);
        CHECK(mat_eq(U * UD, Id(), EXACT));
    }
    SUBCASE("Z") {
        auto U = Z();
        auto UD = adjoint(U);
        CHECK(mat_eq(U * UD, Id(), EXACT));
    }
    SUBCASE("B0") {
        auto U = B0();
        auto UD = adjoint(U);
        INFO("B0 * d\n", U * UD);
        INFO("id\n", Id());
        CHECK_FALSE(mat_eq(U * UD, Id(), EXACT));
    }
    SUBCASE("B1") {
        auto U = B1();
        auto UD = adjoint(U);
        INFO("B1 * d\n", U * UD);
        INFO("id \n", Id());
        CHECK_FALSE(mat_eq(U * UD, Id(), EXACT));
    }
    SUBCASE("CX") {
        auto U = CX();
        auto UD = adjoint(U);
        CHECK(mat_eq(U * UD, Id(2), EXACT));
    }
    SUBCASE("Id(n)") {
        for (int i = 1; i < 10; i++) {
            auto U = Id(i);
            auto UD = adjoint(U);
            CHECK(mat_eq(U * UD, Id(i), EXACT));
        }
    }

    std::vector<double> angles = {0, 1, 45, 90, 180, 270, 360, 720};
    SUBCASE("RX(theta)") {
        for (double ok : angles) {
            auto U = RX(ok);
            auto UD = adjoint(U);
            INFO("RX(", ok, ") is : ", U);
            INFO("RX(", ok, ")† is : ", UD);
            INFO("RX(", ok, ")*RX(", ok, ") is : ", U * UD);
            CHECK(mat_eq(U * UD, Id(), DEC14));
        }
    }
    SUBCASE("RY(theta)") {
        for (double ok : angles) {
            auto U = RY(ok);
            auto UD = adjoint(U);
            INFO("RY(", ok, ") is : ", U);
            INFO("RY(", ok, ")† is : ", UD);
            INFO("RY(", ok, ")*RY(", ok, ") is : ", U * UD);
            CHECK(mat_eq(U * UD, Id(), DEC14));
        }
    }
    SUBCASE("RZ(theta)") {
        for (double ok : angles) {
            auto U = RZ(ok);
            auto UD = adjoint(U);
            INFO("RZ(", ok, ") is : ", U);
            INFO("RZ(", ok, ")† is : ", UD);
            INFO("RZ(", ok, ")*RZ(", ok, ") is : ", U * UD);
            CHECK(mat_eq(U * UD, Id(), DEC14));
        }
    }
    SUBCASE("CX(x,y)") {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if (i != j) {
                    auto U = CG(X(), i, j);
                    auto UD = adjoint(U);
                    auto RU = CG(X(), j, i);
                    auto RUD = adjoint(RU);
                    CHECK(mat_eq(U * UD, Id(slog2(U.n_rows)), DEC14));
                    CHECK(mat_eq(RU * RUD, Id(slog2(U.n_rows)), DEC14));
                }
            }
        }
    }
    SUBCASE("SWAP(x,y)") {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if (i != j) {
                    auto U = SWAP(i, j);
                    auto UD = adjoint(U);
                    auto RU = SWAP(j, i);
                    auto RUD = adjoint(RU);
                    CHECK(mat_eq(U * UD, Id(slog2(U.n_rows)), DEC14));
                    CHECK(mat_eq(RU * RUD, Id(slog2(U.n_rows)), DEC14));
                }
            }
        }
    }
}

TEST_CASE("SWAP Gate Properties") {
    SUBCASE("SWAP is Self-Inverse") {
        for (int q1 = 0; q1 < 3; q1++) {
            for (int q2 = 0; q2 < 3; q2++) {
                if (q1 != q2) {
                    cx_mat swap_gate = SWAP(q1, q2);
                    cx_mat swapped_back = SWAP(q1, q2) * swap_gate;
                    int n = slog2(swapped_back.n_rows);
                    CHECK(mat_eq(swapped_back, Id(n), DEC14));
                }
            }
        }
    }
}

TEST_CASE("slog2 Function") {
    SUBCASE("Power of Two Values") {
        CHECK_EQ(slog2(1), 0);
        CHECK_EQ(slog2(2), 1);
        CHECK_EQ(slog2(4), 2);
        CHECK_EQ(slog2(8), 3);
        CHECK_EQ(slog2(16), 4);
        CHECK_EQ(slog2(32), 5);
    }

    SUBCASE("Non Power of Two Values") {
        // 3 >> 1 = 1, then 1 >> 1 = 0 -> count=1
        CHECK_EQ(slog2(3), 1);
        // 5 >> 1 = 2, then 2 >> 1 = 1, then 1 >> 1 = 0 -> count=2
        CHECK_EQ(slog2(5), 2);
         // 7 >> 1 = 3, then 3 >> 1 = 1, then 1 >> 1 = 0 -> count=2
        CHECK_EQ(slog2(7), 2);
        // 15 >> 1 = 7, then 7 >> 1 = 3, then 3 >> 1 = 1,
        // then 1 >> 1 = 0 -> count=3
        CHECK_EQ(slog2(15), 3);
    }
}
