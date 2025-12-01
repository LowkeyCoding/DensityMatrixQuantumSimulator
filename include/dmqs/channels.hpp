#pragma once
#include <armadillo>

#include <cmath>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <dmqs/gates.hpp>

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

typedef cx_mat::fixed<2, 2> kraus_t;
typedef std::function<vector<kraus_t>(const double&)> channel_t;

vector<kraus_t> amplitude_damping_ops(const double& p);
vector<kraus_t> phase_damping_ops(const double& p);
vector<kraus_t> depolarizing_ops(const double& p);
vector<kraus_t> bit_flip_ops(const double& p);
vector<kraus_t> phase_flip_ops(const double& p);
vector<kraus_t> bit_phase_flip_ops(const double& p);
channel_t u_channel_to_ops_f(u_channel channel);
cx_mat apply_channel(const cx_mat &rho, const vector<kraus_t> &kraus_ops);
vector<kraus_t> generalized_amplitude_damping_ops(const double& p,
                                            const double& gamma);
