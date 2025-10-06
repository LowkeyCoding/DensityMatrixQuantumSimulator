
#include <armadillo>
using namespace arma;
#ifndef GATES_H
#define GATES_H
cx_mat X();
cx_mat CX();
cx_mat Y();
cx_mat Z();
cx_mat H();
cx_mat B0();
cx_mat B1();
cx_mat Id();
cx_mat RX(double theta);
cx_mat RY(double theta);
cx_mat RZ(double theta);
#endif