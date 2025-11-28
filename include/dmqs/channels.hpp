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

typedef cx_mat::fixed<2, 2> qubit;
typedef vector<qubit> kraus_ops;
typedef std::function<kraus_ops(const double&)> channel_f;

kraus_ops amplitude_damping_ops(const double& p);
kraus_ops phase_damping_ops(const double& p);
kraus_ops depolarizing_ops(const double& p);
kraus_ops bit_flip_ops(const double& p);
kraus_ops phase_flip_ops(const double& p);
kraus_ops bit_phase_flip_ops(const double& p);
channel_f u_channel_to_ops_f(u_channel channel);
cx_mat apply_channel(const cx_mat &rho, const kraus_ops &kraus_ops);
kraus_ops generalized_amplitude_damping_ops(const double& p,
                                            const double& gamma);
