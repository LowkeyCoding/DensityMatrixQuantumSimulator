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
cx_mat BinaryStringToDensityMatrix(const string bin);
cx_mat ApplyGateToDensityMatrix(cx_mat matrix, cx_mat gate);
cx_mat ApplyGate(cx_mat rho, u_gate gate, int qubit);
cx_mat ApplyCGate(cx_mat rho, u_gate gate, int target, int control);
cx_mat GateToNQubitSystem(cx_mat gate, int target, int n);
cx_mat ApplyAmplitudeDampeningAndDephasing(cx_mat rho, double* T1, double* T2, double t);
cx_mat PartialTrace(cx_mat rho, vector<int> targets);
int Sample(const cx_mat rho, double random);
int PartialSample(const cx_mat rho, int target, double random);
int PartialSample(const cx_mat rho, vector<int> targets, double random);
cx_mat BasisProjection(const cx_mat rho, int target, int state);
cx_mat BasisProjections(const cx_mat rho, vector<int> targets, int state);
int rearrangeBits(int i, vector<int> a);