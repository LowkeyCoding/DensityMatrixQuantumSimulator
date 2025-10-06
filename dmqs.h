#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <cassert>

using namespace std;
using namespace arma;

cx_mat BinaryStringToDensityMatrix(const string bin);
cx_mat ApplyGateToDensityMatrix(cx_mat matrix, cx_mat gate);
cx_mat GateToNQubitSystem(cx_mat gate, int target, int n);
cx_mat PartialTrace(cx_mat rho, vector<int> targets);
vector<int> Sample(cx_mat rho, vector<double> random_values);
int rearrangeBits(int i, vector<int> a);
bool DensityMatrixApproxEq(cx_mat d1, cx_mat d2, double delta);