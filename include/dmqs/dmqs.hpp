#pragma once
#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <cassert>
#include <gates/gates.hpp>

using namespace std;
using namespace arma;

bool IsPure(cx_mat rho, double delta);
vector<double> GetPropabilities(cx_mat rho);
string DensityMatrixToProbabilityString(cx_mat rho);
cx_mat BinaryStringToDensityMatrix(const string bin);
cx_mat ApplyGateToDensityMatrix(cx_mat matrix, cx_mat gate);
cx_mat ApplyGate(cx_mat rho, u_gate gate, int qubit);
cx_mat ApplyCGate(cx_mat rho, u_gate gate, int target, int control);
cx_mat GateToNQubitSystem(cx_mat gate, int target, int n);
cx_mat ApplyAmplitudeDampeningAndDephasing(cx_mat rho, double* T1, double* T2, double t);
cx_mat PartialTrace(cx_mat rho, vector<int> targets);
int Sample(const cx_mat rho, double random_value);
cx_mat MeasurementGate(const cx_mat rho, int target, double random_value);
int rearrangeBits(int i, vector<int> a);