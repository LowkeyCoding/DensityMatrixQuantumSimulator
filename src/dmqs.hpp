#pragma once
#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <cassert>
#include "gates/gates.hpp"

using namespace std;
using namespace arma;

bool IsPure(cx_mat rho, double delta);
cx_mat BinaryStringToDensityMatrix(const string bin);
cx_mat ApplyGateToDensityMatrix(cx_mat matrix, cx_mat gate);
cx_mat ApplyGate(cx_mat rho, u_gate gate, size_t qubit);
cx_mat GateToNQubitSystem(cx_mat gate, size_t target, size_t n);
cx_mat ApplyAmplitudeDampeningAndDephasing(cx_mat rho, double* T1, double* T2, double t);
cx_mat PartialTrace(cx_mat rho, vector<size_t> targets);
int Sample(const cx_mat rho, double random_value);
size_t rearrangeBits(size_t i, vector<size_t> a);