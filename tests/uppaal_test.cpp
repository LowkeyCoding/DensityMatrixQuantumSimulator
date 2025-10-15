#include "../bindings/uppaal.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../libs/doctest.h"

TEST_CASE("State intializer"){
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
    CHECK(std::equal(std::begin(r00), std::end(r00), std::begin(c00)));
    UInitBinState(r01, qubits, "01");
    CHECK(std::equal(std::begin(r01), std::end(r01), std::begin(c01)));
    UInitBinState(r10, qubits, "10");
    CHECK(std::equal(std::begin(r10), std::end(r10), std::begin(c10)));
    UInitBinState(r11, qubits, "11");
    CHECK(std::equal(std::begin(r11), std::end(r11), std::begin(c11)));

    SUBCASE("sample states"){
        double random[20] = {0.0045, 0.7908, 0.9903, 0.1085, 0.8035, 0.6303, 0.2853, 0.1125, 0.3623, 0.5738, 0.9714, 0.1138, 0.9487, 0.3961, 0.0949, 0.0730, 0.3724, 0.5300, 0.9606, 0.8233};
        int res0[20];
        int res1[20];
        int res2[20];
        int res3[20];
        UMeasureAll(r00,qubits, random, res0, 20);
        for (int i = 0; i < 20; i++) {
            CHECK(res0[i] == 0);
        }
    
        UMeasureAll(r01,qubits, random, res1, 20);
        for (int i = 0; i < 20; i++) {
            CHECK(res1[i] == 1);
        }
        UMeasureAll(r10,qubits, random, res2, 20);
        for (int i = 0; i < 20; i++) {
            CHECK(res2[i] == 2);
        }
        UMeasureAll(r11,qubits, random, res3, 20);
        for (int i = 0; i < 20; i++) {
            CHECK(res3[i] == 3);
        }
    }

}