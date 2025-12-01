#pragma once
#include <armadillo>
#include <cassert>
using std::invalid_argument, std::to_string;

using arma::cx_mat, arma::cx_double, arma::fill::eye;

// 1 qubit gate
typedef cx_mat::fixed<2, 2> gate1_t;

enum u_gate {
  GID,
  GX,
  GY,
  GZ,
  GH,
  GB0,
  GB1,
};

gate1_t X();
cx_mat CX();
gate1_t Y();
gate1_t Z();
gate1_t H();
gate1_t B0();
gate1_t B1();
gate1_t Id();
cx_mat Id(int n);
cx_mat RX(double theta);
cx_mat RY(double theta);
cx_mat RZ(double theta);
cx_mat CG(const cx_mat& gate, int q1, int q2);
cx_mat SWAP(int q1, int q2);

bool mat_eq(const cx_mat& rho1, const cx_mat& rho2, double delta);
cx_mat adjoint(const cx_mat& M);
int slog2(int n);
