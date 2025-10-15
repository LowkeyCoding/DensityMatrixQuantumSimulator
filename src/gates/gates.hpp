#pragma once
#include <armadillo>
#include <cassert>
using namespace arma;

enum u_gate {
  GId,
  GX,
  GY,
  GZ,
  GH,
  GCX,
  GB0,
  GB1,
};

cx_mat X();
cx_mat CX();
cx_mat Y();
cx_mat Z();
cx_mat H();
cx_mat B0();
cx_mat B1();
cx_mat Id();
cx_mat Id(int n);
cx_mat RX(double theta);
cx_mat RY(double theta);
cx_mat RZ(double theta);
cx_mat CG(cx_mat gate, int q1, int q2);

bool DensityMatrixApproxEq(cx_mat rho1, cx_mat rho2, double delta);