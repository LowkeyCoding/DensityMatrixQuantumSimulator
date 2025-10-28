#pragma once
#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <cassert>
#include <gates/gates.hpp>

using namespace std;
using namespace arma;

bool IsPure(const cx_mat& rho, double delta);
cx_mat BinaryStringToDensityMatrix(const string& bin);
cx_mat ApplyGateToDensityMatrix(const cx_mat& rho, const cx_mat& U);
cx_mat ApplyGate(const cx_mat& rho, u_gate gate, int qubit);
cx_mat ApplyCGate(const cx_mat& rho, u_gate gate, int control, int target);
cx_mat GateToNQubitSystem(const cx_mat& U1, int target, int n);
cx_mat ApplyAmplitudeDampeningAndDephasing(const cx_mat& rho, double* T1, double* T2, double t);
cx_mat PartialTrace(const cx_mat& rho, const vector<int>& targets);
int Sample(const cx_mat& rho, double random);
int PartialSample(const cx_mat& rho, int target, double random);
int PartialSample(const cx_mat& rho, const vector<int>& targets, double random);
cx_mat BasisProjection(const cx_mat& rho, int target, int state);
cx_mat BasisProjections(const cx_mat& rho, const vector<int>& targets, int state);
int rearrangeBits(int i, const vector<int>& a);