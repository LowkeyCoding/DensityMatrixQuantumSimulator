#pragma once
#include <armadillo>

#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <cassert>

using namespace arma;
using namespace std;

vector<cx_mat> amplitude_damping_ops(double gamma);
vector<cx_mat> phase_damping_ops(double gamma);
vector<cx_mat> depolarizing_ops(double p);
vector<cx_mat> bit_flip_ops(double p);
vector<cx_mat> phase_flip_ops(double p);
vector<cx_mat> bit_phase_flip_ops(double p);
cx_mat multi_kron(const vector<cx_mat>& matrices);
cx_mat apply_channel(const cx_mat& rho, const vector<cx_mat>& kraus_ops, int n_qubits);