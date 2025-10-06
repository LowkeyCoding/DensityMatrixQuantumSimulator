#include "gates.h"
/// @brief Basis state |0⟩
cx_mat B0(){
     cx_mat::fixed<2, 2> gate = {
        cx_double(1,0),cx_double(0,0),
        cx_double(0,0),cx_double(0,0)
    };
    return gate;
}
/// @brief Basis state |1⟩
cx_mat B1(){
     cx_mat::fixed<2, 2> gate = {
        cx_double(0,0),cx_double(0,0),
        cx_double(0,0),cx_double(1,0),
    };
    return gate;
}
/// @brief Identity gate
cx_mat Id() {
    cx_mat::fixed<2, 2> gate = {
        cx_double(1,0),cx_double(0,0),
        cx_double(0,0),cx_double(1,0),
    };
    return gate;
}

/// @brief X gate
cx_mat X() {
    cx_mat::fixed<2, 2> gate = {
        cx_double(0,0),cx_double(1,0),
        cx_double(1,0),cx_double(0,0),
    };
    return gate;
} 

/// @brief Y gate 
cx_mat Y() {
    cx_mat::fixed<2, 2> gate = {
        cx_double(0,0),cx_double(0,-1),
        cx_double(0,1),cx_double(0,0),
    };
    return gate;
}

/// @brief Z gate 
cx_mat Z() {
    cx_mat::fixed<2, 2> gate = {
        cx_double(1,0),cx_double(0,0),
        cx_double(0,0),cx_double(-1,0),
    };
    return gate;
}

/// @brief Hadamard gate 
cx_mat H() {
    cx_mat::fixed<2, 2> gate = {
        cx_double(1/sqrt(2),0),cx_double(1/sqrt(2),0),
        cx_double(1/sqrt(2),0),cx_double(-(1/sqrt(2)),0),
    };
    return gate;
}
/// @brief Controlled X gate 
cx_mat CX() {
    cx_mat::fixed<4, 4> gate = {
        cx_double(1,0),cx_double(0,0), cx_double(0,0),cx_double(0,0),
        cx_double(0,0),cx_double(1,0), cx_double(0,0),cx_double(0,0),
        cx_double(0,0),cx_double(0,0), cx_double(0,0),cx_double(1,0),
        cx_double(0,0),cx_double(0,0), cx_double(1,0),cx_double(0,0),
    };
    return gate;
}
/// @brief RX gate (Rotates the bloch sphere theta degrees around the X axis)
/// @param theta Angle in degrees
/// @return RX gate with a rotation theta
cx_mat RX(double theta){
    double theta_r = (theta / 180.0 * M_PI) / 2; // Half the radians of theta
    cx_mat::fixed<2,2> rx_theta = {
        cx_double(cos(theta_r/2),0), cx_double(0,-sin(theta_r)),
        cx_double(0,-sin(theta_r/2)), cx_double(cos(theta_r/2),0)
    };
    return rx_theta.reshape(2,2);
}


/// @brief RY gate (Rotates the bloch sphere theta degrees around the Y axis)
/// @param theta Angle in degrees
/// @return RY gate with a rotation theta
cx_mat RY(double theta){
    double theta_r = (theta / 180.0 * M_PI) / 2; // Half the radians of theta
    cx_mat::fixed<2,2> ry_theta = {
        cx_double(cos(theta_r/2),0), cx_double(-sin(theta_r),0),
        cx_double(sin(theta_r/2),0), cx_double(cos(theta_r/2),0)
    };
    return ry_theta.reshape(2,2);
}

/// @brief RZ gate (Rotates the bloch sphere theta degrees around the Z axis)
/// @param theta Angle in degrees
/// @return RZ gate with a rotation theta
cx_mat RZ(double theta){
    double theta_r = (theta / 180.0 * M_PI) / 2; // Half the radians of theta
    cx_mat::fixed<2,2> rz_theta = {
        cx_double(cos(theta_r/2),-sin(theta_r/2)), cx_double(0,0),
        cx_double(0,0), cx_double(cos(theta_r/2),sin(theta_r/2))
    };
    return rz_theta;
}