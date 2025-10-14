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

cx_mat Id(int n) {
    return cx_mat(pow(2,n),pow(2,n),fill::eye);
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

/// @brief Takes an abitrary 1 qubit gate and control and target qubits and constructs a controlled version of the gate.
/// @param gate The gate to apply to the target qubit when the control is 1.
/// @param control The qubit that controls whether the gate should be applied.
/// @param target
/// @return Controlled gate
cx_mat CG(cx_mat gate, int control, int target) {
    cx_mat mat_control = B0();
    cx_mat mat_target = B1();
    int dist_between_qubits = abs(control-target);

    if (dist_between_qubits > 1) {
        if (control > target) {
            mat_target = kron(Id(dist_between_qubits-1), mat_target);
        } else {
            mat_target = kron(mat_target, Id(dist_between_qubits-1));
        }
    }

    if(control > target) {
       mat_control = kron(Id(dist_between_qubits), mat_control);
       mat_target =  kron(gate, mat_target);
    } else {
       mat_target = kron(mat_target, gate);
       mat_control = kron(mat_control,Id(dist_between_qubits));
    }

    return mat_control + mat_target;
}