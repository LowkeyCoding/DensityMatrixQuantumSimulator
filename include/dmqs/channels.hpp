#pragma once
#include <armadillo>

#include <cmath>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include "gates/gates.hpp"

using std::vector, std::begin, std::end;
using arma::cx_mat, arma::cx_double, arma::fill::zeros, arma::fill::ones;

enum u_channel {
    AMPLITUDE_DAMPING,
    PHASE_DAMPING,
    BIT_FLIP,
    PHASE_FLIP,
    BIT_PHASE_FLIP,
    DEPOLARIZING,
};

vector<cx_mat> amplitude_damping_ops(double p);
vector<cx_mat> phase_damping_ops(double p);
vector<cx_mat> depolarizing_ops(double p);
vector<cx_mat> bit_flip_ops(double p);
vector<cx_mat> phase_flip_ops(double p);
vector<cx_mat> bit_phase_flip_ops(double p);
cx_mat multi_kron(const vector<cx_mat> &matrices);
std::function<vector<cx_mat>(double)> u_channel_to_ops_f(u_channel channel);
cx_mat apply_channel(const cx_mat &rho, const vector<cx_mat> &kraus_ops);
cx_mat apply_channel_with_ops_for_each_qubit(const cx_mat &rho,
                                             const vector<cx_mat> &kraus_ops);
