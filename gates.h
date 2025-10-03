
#include <armadillo>
using namespace arma;
#ifndef GATES_H
#define GATES_H
/// @brief Basis state |0⟩
cx_mat::fixed<2, 2> S0 = {
    cx_double(1,0),cx_double(0,0),
    cx_double(0,0),cx_double(0,0)
};
/// @brief Basis state |1⟩
cx_mat::fixed<2, 2> S1 = {
    cx_double(0,0),cx_double(0,0),
    cx_double(0,0),cx_double(1,0),
};
/// @brief Identity gate
cx_mat::fixed<2, 2> Id = {
    cx_double(1,0),cx_double(0,0),
    cx_double(0,0),cx_double(1,0),
};

/// @brief X gate 
cx_mat::fixed<2, 2> X = {
    cx_double(0,0),cx_double(1,0),
    cx_double(1,0),cx_double(0,0),
};

/// @brief Y gate 
cx_mat::fixed<2, 2> Y = {
    cx_double(0,0),cx_double(0,-1),
    cx_double(0,1),cx_double(0,0),
};
/// @brief Z gate 
cx_mat::fixed<2, 2> Z = {
    cx_double(1,0),cx_double(0,0),
    cx_double(0,0),cx_double(-1,0),
};
/// @brief Hadamard gate 
cx_mat::fixed<2, 2> H = {
    cx_double(1/sqrt(2),0),cx_double(1/sqrt(2),0),
    cx_double(1/sqrt(2),0),cx_double(-(1/sqrt(2)),0),
};
/// @brief Controlled X gate 
cx_mat::fixed<4, 4> CX = {
    cx_double(1,0),cx_double(0,0), cx_double(0,0),cx_double(0,0),
    cx_double(0,0),cx_double(1,0), cx_double(0,0),cx_double(0,0),
    cx_double(0,0),cx_double(0,0), cx_double(0,0),cx_double(1,0),
    cx_double(0,0),cx_double(0,0), cx_double(1,0),cx_double(0,0),
};
/// @brief RX gate (Rotates the bloch sphere theta degrees around the X axis)
/// @param theta Angle in degrees
/// @return RX gate with a rotation theta
static cx_mat RX(double theta){
    double theta_r = (theta / 180.0 * M_PI) / 2; // Half the radians of theta
    cx_mat rx_theta = {
        cx_double(cos(theta_r/2),0), cx_double(0,-sin(theta_r)),
        cx_double(0,-sin(theta_r/2)), cx_double(cos(theta_r/2),0)
    };
    return rx_theta.reshape(2,2);
}


/// @brief RY gate (Rotates the bloch sphere theta degrees around the Y axis)
/// @param theta Angle in degrees
/// @return RY gate with a rotation theta
static cx_mat RY(double theta){
    double theta_r = (theta / 180.0 * M_PI) / 2; // Half the radians of theta
    cx_mat ry_theta = {
        cx_double(cos(theta_r/2),0), cx_double(-sin(theta_r),0),
        cx_double(sin(theta_r/2),0), cx_double(cos(theta_r/2),0)
    };
    return ry_theta.reshape(2,2);
}

/// @brief RZ gate (Rotates the bloch sphere theta degrees around the Z axis)
/// @param theta Angle in degrees
/// @return RZ gate with a rotation theta
static cx_mat RZ(double theta){
    double theta_r = (theta / 180.0 * M_PI) / 2; // Half the radians of theta
    cx_mat rz_theta = {
        cx_double(cos(theta_r/2),-sin(theta_r/2)), cx_double(0,0),
        cx_double(0,0), cx_double(cos(theta_r/2),sin(theta_r/2))
    };
    return rz_theta.reshape(2,2);
}

/// @brief RV gate (Rotats all axis theta amount of degrees)
/// @param theta_x 
/// @param theta_y 
/// @param theta_z 
/// @return 
static cx_mat RV(double theta_x, double theta_y, double theta_z) {
    double magnitue_sqr = theta_x*theta_x + theta_y*theta_y + theta_z*theta_z;
    if (magnitue_sqr == 0) {
        return Id;
    }
    double radians = sqrt(magnitue_sqr);
    theta_x /= radians;
    theta_y /= radians;
    theta_z /= radians;
    double sine = sin(radians/2);
    double cosine = cos(radians/2);
    cx_mat rv_rot = {
        cx_double(cosine, -theta_z*sine),       cx_double(-theta_y*sine,-theta_x*sine),
        cx_double(theta_y*sine,-theta_x*sine),  cx_double(cosine, theta_z*sine)
    };
    return rv_rot.reshape(2,2);
}
#endif