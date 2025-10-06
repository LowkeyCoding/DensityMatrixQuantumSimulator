#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <cassert>
#include "gates.h"

using namespace std;
using namespace arma;

cx_mat BinaryStringToDensityMatrix(const string bin);
cx_mat ApplyGateToDensityMatrix(cx_mat matrix, cx_mat gate);
cx_mat GateToNQubitSystem(cx_mat gate, size_t target, size_t n);
cx_mat ApplyGate(cx_mat rho, string US, size_t qubit);
cx_mat PartialTrace(cx_mat rho, vector<size_t> targets);
vector<size_t> Sample(cx_mat rho, vector<double> random_values);
size_t rearrangeBits(size_t i, vector<size_t> a);
bool DensityMatrixApproxEq(cx_mat d1, cx_mat d2, double delta);