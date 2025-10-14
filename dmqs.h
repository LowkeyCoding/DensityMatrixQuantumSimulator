#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <cassert>
#include "gates.h"

using namespace std;
using namespace arma;
#ifndef DMQS_H
#define DMQS_H

bool IsPure(cx_mat rho, double delta);
cx_mat BinaryStringToDensityMatrix(const string bin);
cx_mat ApplyGateToDensityMatrix(cx_mat matrix, cx_mat gate);
cx_mat ApplyGate(cx_mat rho, u_gate gate, size_t qubit);
cx_mat GateToNQubitSystem(cx_mat gate, size_t target, size_t n);
cx_mat ApplyAmplitudeDampeningAndDephasing(cx_mat rho, double* T1, double* T2, double t);
cx_mat PartialTrace(cx_mat rho, vector<size_t> targets);
vector<int32_t> Sample(cx_mat rho, vector<double> random_values);
size_t rearrangeBits(size_t i, vector<size_t> a);
bool DensityMatrixApproxEq(cx_mat d1, cx_mat d2, double delta);

#endif