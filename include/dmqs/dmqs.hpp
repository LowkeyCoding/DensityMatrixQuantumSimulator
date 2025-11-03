#pragma once
#include <armadillo>

#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <cassert>
#include <gates/gates.hpp>

using std::vector, std::string, std::invalid_argument, std::to_string,
      std::max, std::min, std::greater, std::map, std::equal,
      std::end, std::begin;

using arma::cx_mat, arma::cx_double, arma::cumsum, arma::any, arma::vec,
      arma::uword, arma::fill::zeros;

bool IsPure(const cx_mat& rho, double delta);
cx_mat BinaryStringToDensityMatrix(const string& bin);
cx_mat ApplyGateToDensityMatrix(const cx_mat& rho, const cx_mat& U);
cx_mat ApplyGate(const cx_mat& rho, u_gate gate, int qubit);
cx_mat ApplyCGate(const cx_mat& rho, u_gate gate, int control, int target);
cx_mat GateToNQubitSystem(const cx_mat& U1, int target, int n);
cx_mat ApplyAmplitudeDampeningAndDephasing(const cx_mat& rho, const double* T1,
                                           const double* T2, double t);
cx_mat PartialTrace(const cx_mat& rho, const vector<int>& targets);
int Sample(const cx_mat& rho, double random);
int PartialSample(const cx_mat& rho, int target, double random);
int PartialSample(const cx_mat& rho, const vector<int>& targets, double random);
cx_mat BasisProjection(const cx_mat& rho, int target, int state);
cx_mat BasisProjections(const cx_mat& rho, const vector<int>& targets,
                        int state);
int rearrangeBits(int i, const vector<int>& a);
