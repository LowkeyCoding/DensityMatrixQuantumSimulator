#include <doctest/doctest.h>
#include <vector>
#include <string>
#include <dmqs/dmqs.hpp>
#include <dmqs/channels.hpp>

#define EXACT 0.0
#define DEC14 1e-14

TEST_CASE("Amplitude dampning") {
    const cx_mat rho_00 = BinaryStringToDensityMatrix("00");
    const cx_mat rho_01 = BinaryStringToDensityMatrix("01");
    const cx_mat rho_10 = BinaryStringToDensityMatrix("10");
    const cx_mat rho_11 = BinaryStringToDensityMatrix("11");
    const cx_mat rho_pp = BinaryStringToDensityMatrix("++");
    const cx_mat rho_mm = BinaryStringToDensityMatrix("--");
    SUBCASE("gamma=0.0") {
        double gamma = 0.0;
        kraus_ops channel = amplitude_damping_ops(1 - gamma);
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_01, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_10);
            CHECK(mat_eq(res, rho_10, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_11);
            CHECK(mat_eq(res, rho_11, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_pp);
            CHECK(mat_eq(res, rho_pp, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_mm);
            CHECK(mat_eq(res, rho_mm, DEC14));
        }
    }
    SUBCASE("gamma=0.25") {
        double gamma = 0.25;
        kraus_ops channel = amplitude_damping_ops(1 - gamma);
        cx_mat expected = cx_mat(4, 4, zeros);
        SUBCASE("Rho=|00><00|") {
            // Should remain unchanged as when gamma goes towards 1,
            // limit gamma->1 means no damping
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            expected(0, 0) = 0.25; // 1 - gamma
            expected(1, 1) = 0.75; // 1 - (1 - gamma)
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            expected(0, 0) = 0.25; // 1 - gamma
            expected(2, 2) = 0.75; // 1 - (1 - gamma)
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            // pow(0.25, 2)
            expected(0, 0) = 0.0625;
            // sqrt(0.75)*sqrt(0.25)*sqrt(0.75)*sqrt(0.25)
            expected(1, 1) = 0.1875;
            // sqrt(0.75)*sqrt(0.25)*sqrt(0.75)*sqrt(0.25)
            expected(2, 2) = 0.1875;
            // pow(0.75, 2)
            expected(3, 3) = 0.5625;
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        /*SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            expected(0, 0) = 25.0/64.0;
            expected(1, 1) = 0.25 * sqrt(0.75) * 0.25 * sqrt(75) * 0.5;
            expected(2, 2) = 0.25 * sqrt(0.75) * 0.25 * sqrt(75) * 0.5;
            expected(3, 3) = pow(sqrt(0.75), 4) * 0.25;
            double a = 0.1875;//sqrt(0.75)*sqrt(0.25)*sqrt(0.75)*sqrt(0.25);
            expected(3, 0) = a;
            expected(2, 1) = a;
            expected(1, 2) = a;
            expected(0, 3) = a;
            double b = 0.2706329386826369; // sqrt(0.75)*0.25*0.25
            expected(1, 0) = b;
            expected(2, 0) = b;
            expected(0, 1) = b;
            expected(0, 2) = b;
            double c = sqrt(0.75)*0.25*0.75; // sqrt(0.75)*0.25*0.75
            expected(3, 1) = c;
            expected(1, 3) = c;
            expected(3, 2) = c;
            expected(2, 3) = c;

            INFO("Rho: \n", rho_pp);
            INFO("Result:\n", res, "\nExpected:\n", expected);
            INFO("Diff:\n", res - expected);
            INFO(std::setprecision(16), res(0, 1));
            CHECK(mat_eq(res, rho_00, DEC14));
        }*/
    }
    SUBCASE("gamma=0.5") {
        double gamma = 0.5;
        kraus_ops channel = amplitude_damping_ops(1 - gamma);
        cx_mat res = apply_channel(rho_01, channel);
        cx_mat expected = cx_mat(4, 4, zeros);
        expected(0, 0) = 0.5;
        expected(1, 1) = 0.5;
        INFO("Result:\n", res, "\nExpected:\n", expected);
        CHECK(mat_eq(res, expected, DEC14));
    }
    SUBCASE("gamma=0.75") {
        double gamma = 0.75;
        kraus_ops channel = amplitude_damping_ops(1 - gamma);
        cx_mat res = apply_channel(rho_01, channel);
        cx_mat expected = cx_mat(4, 4, zeros);
        expected(0, 0) = 0.75;
        expected(1, 1) = 0.25;
        INFO("Result:\n", res, "\nExpected:\n", expected);
        CHECK(mat_eq(res, expected, DEC14));
    }
    SUBCASE("gamma=1.0") {
        double gamma = 1.0;
        kraus_ops channel = amplitude_damping_ops(1 - gamma);
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
    }
}

TEST_CASE("Phase Dampning") {
    const cx_mat rho_00 = BinaryStringToDensityMatrix("00");
    const cx_mat rho_01 = BinaryStringToDensityMatrix("01");
    const cx_mat rho_10 = BinaryStringToDensityMatrix("10");
    const cx_mat rho_11 = BinaryStringToDensityMatrix("11");
    const cx_mat rho_pp = BinaryStringToDensityMatrix("++");
    const cx_mat rho_mm = BinaryStringToDensityMatrix("--");

    SUBCASE("gamma=0.0") {
        double gamma = 0.0;
        kraus_ops channel = phase_damping_ops(1 - gamma);
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_01);
            CHECK(mat_eq(res, rho_01, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_10);
            CHECK(mat_eq(res, rho_10, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_11);
            CHECK(mat_eq(res, rho_11, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_pp);
            CHECK(mat_eq(res, rho_pp, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_mm);
            CHECK(mat_eq(res, rho_mm, DEC14));
        }
    }

    SUBCASE("gamma=0.25") {
        double gamma = 0.25;
        kraus_ops channel = phase_damping_ops(1 - gamma);
        INFO("Channel Kraus Operators:\n");
        for (const cx_mat& K : channel) {
            INFO(K, "\n");
        }
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_01);
            CHECK(mat_eq(res, rho_01, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_10);
            CHECK(mat_eq(res, rho_10, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_11);
            CHECK(mat_eq(res, rho_11, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;
            // 0.140625 = pow(0.75, 2) * 0.25 or pow(sqrt(0.75), 4) * 0.25
            double a = pow(0.75, 2) * 0.25;
            expected(3, 0) = a;
            expected(2, 1) = a;
            expected(1, 2) = a;
            expected(0, 3) = a;
            // 0.1875 = sqrt(0.75) * sqrt(0.25) * sqrt(0.75) * sqrt(0.25)
            double b = sqrt(0.75) * sqrt(0.25) * sqrt(0.75) * sqrt(0.25);
            expected(0, 1) = b;
            expected(0, 2) = b;
            expected(1, 0) = b;
            expected(1, 3) = b;
            expected(2, 0) = b;
            expected(2, 3) = b;
            expected(3, 1) = b;
            expected(3, 2) = b;
            INFO("Rho: \n", rho_pp);
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;
            // 0.140625 = pow(0.75, 2) * 0.25 or pow(sqrt(0.75), 4) * 0.25
            double a = pow(0.75, 2) * 0.25;
            expected(3, 0) = a;
            expected(2, 1) = a;
            expected(1, 2) = a;
            expected(0, 3) = a;
            // -0.1875 = -sqrt(0.75) * sqrt(0.25) * sqrt(0.75) * sqrt(0.25)
            double b = -sqrt(0.75) * sqrt(0.25) * sqrt(0.75) * sqrt(0.25);
            expected(0, 1) = b;
            expected(0, 2) = b;
            expected(1, 0) = b;
            expected(1, 3) = b;
            expected(2, 0) = b;
            expected(2, 3) = b;
            expected(3, 1) = b;
            expected(3, 2) = b;
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
    }

    SUBCASE("gamma=0.5") {
        double gamma = 0.5;
        kraus_ops channel = phase_damping_ops(1 - gamma);
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_01);
            CHECK(mat_eq(res, rho_01, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_10);
            CHECK(mat_eq(res, rho_10, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_11);
            CHECK(mat_eq(res, rho_11, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;

            expected(3, 0) = 0.25 * 0.25;
            expected(2, 1) = 0.25 * 0.25;
            expected(1, 2) = 0.25 * 0.25;
            expected(0, 3) = 0.25 * 0.25;

            expected(0, 1) = 0.25 * 0.5;
            expected(0, 2) = 0.25 * 0.5;
            expected(1, 0) = 0.25 * 0.5;
            expected(1, 3) = 0.25 * 0.5;
            expected(2, 0) = 0.25 * 0.5;
            expected(2, 3) = 0.25 * 0.5;
            expected(3, 1) = 0.25 * 0.5;
            expected(3, 2) = 0.25 * 0.5;

            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;

            expected(3, 0) = 0.25 * 0.25;
            expected(2, 1) = 0.25 * 0.25;
            expected(1, 2) = 0.25 * 0.25;
            expected(0, 3) = 0.25 * 0.25;

            expected(0, 1) = -0.25 * 0.5;
            expected(0, 2) = -0.25 * 0.5;
            expected(1, 0) = -0.25 * 0.5;
            expected(1, 3) = -0.25 * 0.5;
            expected(2, 0) = -0.25 * 0.5;
            expected(2, 3) = -0.25 * 0.5;
            expected(3, 1) = -0.25 * 0.5;
            expected(3, 2) = -0.25 * 0.5;

            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
    }

    SUBCASE("gamma=0.75") {
        double gamma = 0.75;
        kraus_ops channel = phase_damping_ops(1 - gamma);
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_01);
            CHECK(mat_eq(res, rho_01, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_10);
            CHECK(mat_eq(res, rho_10, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_11);
            CHECK(mat_eq(res, rho_11, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;
            double a = pow(0.25, 3);
            expected(3, 0) = a;
            expected(2, 1) = a;
            expected(1, 2) = a;
            expected(0, 3) = a;
            double b = pow(0.25, 2);
            expected(0, 1) = b;
            expected(0, 2) = b;
            expected(1, 0) = b;
            expected(1, 3) = b;
            expected(2, 0) = b;
            expected(2, 3) = b;
            expected(3, 1) = b;
            expected(3, 2) = b;
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;
            double a = pow(0.25, 3);
            expected(3, 0) = a;
            expected(2, 1) = a;
            expected(1, 2) = a;
            expected(0, 3) = a;
            double b = -pow(0.25, 2);
            expected(0, 1) = b;
            expected(0, 2) = b;
            expected(1, 0) = b;
            expected(1, 3) = b;
            expected(2, 0) = b;
            expected(2, 3) = b;
            expected(3, 1) = b;
            expected(3, 2) = b;
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
    }

    SUBCASE("gamma=1.0") {
        double gamma = 1.0;
        kraus_ops channel = phase_damping_ops(1 - gamma);
        SUBCASE("Rho=|00><00|") {
            cx_mat res = apply_channel(rho_00, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_00);
            CHECK(mat_eq(res, rho_00, DEC14));
        }
        SUBCASE("Rho=|01><01|") {
            cx_mat res = apply_channel(rho_01, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_01);
            CHECK(mat_eq(res, rho_01, DEC14));
        }
        SUBCASE("Rho=|10><10|") {
            cx_mat res = apply_channel(rho_10, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_10);
            CHECK(mat_eq(res, rho_10, DEC14));
        }
        SUBCASE("Rho=|11><11|") {
            cx_mat res = apply_channel(rho_11, channel);
            INFO("Result:\n", res, "\nExpected:\n", rho_11);
            CHECK(mat_eq(res, rho_11, DEC14));
        }
        SUBCASE("Rho=|++><++|") {
            cx_mat res = apply_channel(rho_pp, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
        SUBCASE("Rho=|--><--|") {
            cx_mat res = apply_channel(rho_mm, channel);
            cx_mat expected = cx_mat(4, 4, zeros);
            expected(0, 0) = 0.25;
            expected(1, 1) = 0.25;
            expected(2, 2) = 0.25;
            expected(3, 3) = 0.25;
            INFO("Result:\n", res, "\nExpected:\n", expected);
            CHECK(mat_eq(res, expected, DEC14));
        }
    }
}
