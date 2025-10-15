#include "../bindings/uppaal.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../libs/doctest.h"

TEST_CASE("State intializer"){
    int qubits = 2;
    double res[32] = {0};
    cx_double d0 = cx_double(0,0);
    cx_double d1 = cx_double(1,0);
    
    double c00[32] = {0};
    c00[0] = 1;
    double c01[32] = {0};
    c01[10] = 1;
    double c10[32] = {0};
    c10[20] = 1;
    double c11[32] = {0};
    c11[30] = 1;
    
    UInitBinState(res, qubits, "00\0");
    CHECK(std::equal(std::begin(res), std::end(res), std::begin(c00)));
}