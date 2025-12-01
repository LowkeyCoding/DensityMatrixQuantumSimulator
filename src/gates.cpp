#include <dmqs/gates.hpp>

/// @brief Basis state |0⟩
gate1_t B0() {
    static const gate1_t gate = {
        {cx_double(1, 0), cx_double(0, 0)},
        {cx_double(0, 0), cx_double(0, 0)}
    };
    return gate;
}

/// @brief Basis state |1⟩
gate1_t B1() {
    static const gate1_t gate = {
        {cx_double(0, 0), cx_double(0, 0)},
        {cx_double(0, 0), cx_double(1, 0)},
    };
    return gate;
}

/// @brief Identity gate
gate1_t Id() {
    static const gate1_t gate = {
        {cx_double(1, 0), cx_double(0, 0)},
        {cx_double(0, 0), cx_double(1, 0)},
    };
    return gate;
}

cx_mat Id(int n) {
    if (n < 1) {
        throw invalid_argument(
            "N must be at least 1 when generating a Id gate");
    }
    return cx_mat(1 << n, 1 << n, eye);
}

/// @brief X gate
gate1_t X() {
    static const gate1_t gate = {
        {cx_double(0, 0), cx_double(1, 0)},
        {cx_double(1, 0), cx_double(0, 0)},
    };
    return gate;
}

/// @brief Y gate
gate1_t Y() {
    static const gate1_t gate = {
        {cx_double(0, 0), cx_double(0, -1)},
        {cx_double(0, 1), cx_double(0, 0)},
    };
    return gate;
}

/// @brief Z gate
gate1_t Z() {
    static const gate1_t gate = {
        {cx_double(1, 0), cx_double(0, 0)},
        {cx_double(0, 0), cx_double(-1, 0)},
    };
    return gate;
}

/// @brief Hadamard gate
gate1_t H() {
    static const gate1_t gate = {
        {cx_double(1/sqrt(2), 0), cx_double(1/sqrt(2), 0)},
        {cx_double(1/sqrt(2), 0), cx_double(-(1/sqrt(2)), 0)},
    };
    return gate;
}

/// @brief Controlled X gate
cx_mat CX() {
    static const cx_mat::fixed<4, 4> gate = {
        {cx_double(1, 0), cx_double(0, 0), cx_double(0, 0), cx_double(0, 0)},
        {cx_double(0, 0), cx_double(1, 0), cx_double(0, 0), cx_double(0, 0)},
        {cx_double(0, 0), cx_double(0, 0), cx_double(0, 0), cx_double(1, 0)},
        {cx_double(0, 0), cx_double(0, 0), cx_double(1, 0), cx_double(0, 0)},
    };
    return gate;
}

/// @brief RX gate (Rotates the bloch sphere theta degrees around the X axis)
/// @param theta Angle in degrees
/// @return RX gate with a rotation theta
cx_mat RX(double theta) {
    double theta_r = theta * M_PI / 180.0 / 2.0; // Convert to half radians

    cx_mat::fixed<2, 2> rx_theta = {
        {cx_double(cos(theta_r), 0), cx_double(0, -sin(theta_r))},
        {cx_double(0, -sin(theta_r)), cx_double(cos(theta_r), 0)}
    };
    return rx_theta;
}

/// @brief RY gate (Rotates the bloch sphere theta degrees around the Y axis)
/// @param theta Angle in degrees
/// @return RY gate with a rotation theta
cx_mat RY(double theta) {
    double theta_r = (theta / 180.0 * M_PI) / 2.0; // Convert to half radians
    cx_mat::fixed<2, 2> ry_theta = {
        {cx_double(cos(theta_r), 0), cx_double(-sin(theta_r), 0)},
        {cx_double(sin(theta_r), 0), cx_double(cos(theta_r), 0)}
    };
    return ry_theta;
}

/// @brief RZ gate (Rotates the bloch sphere theta degrees around the Z axis)
/// @param theta Angle in degrees
/// @return RZ gate with a rotation theta
cx_mat RZ(double theta) {
    double theta_r = (theta / 180.0 * M_PI) / 2.0; // Convert to half radians
    cx_mat::fixed<2, 2> rz_theta = {
        {std::exp(cx_double(0, -theta_r)), cx_double(0, 0)},
        {cx_double(0, 0), std::exp(cx_double(0, theta_r))},
    };
    return rz_theta;
}

/// @brief Takes an abitrary 1 qubit gate and control and target qubits
///        and constructs a controlled version of the gate.
/// @param gate The gate to apply to the target qubit when the control is 1.
/// @param control The qubit that controls whether the gate should be applied.
/// @param target
/// @return Controlled gate
cx_mat CG(const cx_mat& gate, int control, int target) {
    if (control == target) {
        throw invalid_argument("Control and target qubit must be different");
    }
    cx_mat mat_control = B0();
    cx_mat mat_target = B1();
    int dist_between_qubits = abs(control - target);

    if (dist_between_qubits > 1) {
        if (control > target) {
            mat_target = kron(Id(dist_between_qubits - 1), mat_target);
        } else {
            mat_target = kron(mat_target, Id(dist_between_qubits - 1));
        }
    }

    if (control > target) {
       mat_control = kron(Id(dist_between_qubits), mat_control);
       mat_target = kron(gate, mat_target);
    } else {
       mat_target = kron(mat_target, gate);
       mat_control = kron(mat_control, Id(dist_between_qubits));
    }

    return mat_control + mat_target;
}

/// @brief SWAP Gate
/// @param q1 The index of the first qubit
/// @param q2 The index of the second qubit
/// @return A swap gate that swaps the value of q1 abd q2.
cx_mat SWAP(int q1, int q2) {
    if (q1 == q2) {
        throw invalid_argument("The qubits to swap has to be different");
    }
    return CG(X(), q1, q2) * CG(X(), q2, q1) * CG(X(), q1, q2);
}

/// @brief Hermitian adjoint (dagger)
/// @param M The input matrix
/// @return M^†
cx_mat adjoint(const cx_mat& M) {
    return conj(M).st();
}

/// @brief Checks equality between two density matrixes
/// @param rho1 First density matrix
/// @param rho2 Second density matrix
/// @param delta Maximum difference allowed for equality
/// @return Whether the two density matrixes are equal
bool mat_eq(const cx_mat& rho1, const cx_mat& rho2, double delta) {
    return approx_equal(rho1, rho2, "absdiff", delta);
}

/// @brief Gets log2 of a number
/// @param n
/// @return log2 of number without decimal part
int slog2(int n) {
    int result = 0;
    while (n >>= 1) result++;
    return result;
}
