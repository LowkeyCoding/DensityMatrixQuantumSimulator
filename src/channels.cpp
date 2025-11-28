#include <dmqs/channels.hpp>
#include <vector>
// Amplitude damping channel Kraus operators
kraus_ops amplitude_damping_ops(const double& p) {
    kraus_ops ops = {
        { {cx_double(1, 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(p), 0)} },
        { {cx_double(0, 0), cx_double(sqrt(1 - p), 0)},
          {cx_double(0, 0), cx_double(0, 0)}}
    };
    return ops;
}

kraus_ops generalized_amplitude_damping_ops(const double& p,
                                            const double& gamma) {
    return {
        { {cx_double(sqrt(p), 0), 0},
          {0, cx_double(sqrt(p)*sqrt(gamma), 0) } },
        { {0, cx_double(sqrt(p)*sqrt(gamma), 0)},
          {0, 0} },
        { {cx_double(sqrt(1-p)*sqrt(1-gamma), 0), 0 },
          {0, cx_double(sqrt(1-p), 0) } },
        { {0, 0},
          {cx_double(sqrt(1-p)*sqrt(gamma), 0), 0} }
    };
}

// Phase damping channel Kraus operators
kraus_ops phase_damping_ops(const double& p) {
    kraus_ops ops = {
        { {cx_double(sqrt(p), 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(p), 0)} },
        { {cx_double(sqrt(1 - p), 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(0, 0)} },
        { {cx_double(0, 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(1 - p), 0)} }
    };
    return ops;
}

// Bit flip channel Kraus operators
kraus_ops bit_flip_ops(const double& p) {
    kraus_ops ops;
    auto E0 = qubit(zeros);
    E0(0, 1) = 1 * sqrt(1 - p);
    E0(1, 0) = 1 * sqrt(1 - p);
    auto E1 = qubit(zeros);
    E1(0, 0) = sqrt(p);
    E1(1, 1) = sqrt(p);
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

kraus_ops phase_flip_ops(const double& p) {
    kraus_ops ops;
    auto E0 = qubit(zeros);
    E0(0, 0) = sqrt(p);
    E0(1, 1) = sqrt(p);
    auto E1 = qubit(zeros);
    E1(0, 0) = sqrt(1 - p);
    E1(1, 1) = sqrt(1 - p);
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

// Bit Phase flip channel Kraus operators
kraus_ops bit_phase_flip_ops(const double& p) {
    kraus_ops ops;
    auto E0 = qubit(zeros);
    E0(0, 0) = sqrt(p);
    E0(1, 1) = sqrt(p);
    auto E1 = qubit(zeros);
    E1(0, 1) = cx_double(0, -sqrt(1 - p));
    E1(1, 0) = cx_double(0, sqrt(1 - p));
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

kraus_ops depolarizing_ops(const double& p) {
    kraus_ops ops;
    double sqrt_p_over_3 = sqrt(p / 3.0);
    qubit E0 = { {cx_double(sqrt(1-p), 0), cx_double(0, 0)},
                 {cx_double(0, 0), cx_double(sqrt(1-p), 0)} };
    qubit E1 = { {cx_double(0, 0), cx_double(0, sqrt_p_over_3)},
                 {cx_double(0, sqrt_p_over_3), cx_double(0, 0)} };
    qubit E2 = { {cx_double(0, 0), cx_double(0, -sqrt_p_over_3)},
                 {cx_double(0, sqrt_p_over_3), cx_double(0, 0)} };
    qubit E3 = { {cx_double(sqrt_p_over_3, 0), cx_double(0, 0)},
                 {cx_double(0, 0), cx_double(-sqrt_p_over_3, 0)} };
    ops.push_back(E0);
    ops.push_back(E1);
    ops.push_back(E2);
    ops.push_back(E3);
    return ops;
}

channel_f u_channel_to_ops_f(u_channel channel) {
    switch (channel) {
        case AMPLITUDE_DAMPING:
            return amplitude_damping_ops;
        case PHASE_DAMPING:
            return phase_damping_ops;
        case BIT_FLIP:
            return bit_flip_ops;
        case PHASE_FLIP:
            return phase_flip_ops;
        case BIT_PHASE_FLIP:
            return bit_phase_flip_ops;
        case DEPOLARIZING:
            return depolarizing_ops;
        default:
            throw std::invalid_argument("Unknown channel type");
    }
}

cx_mat apply_channel(const cx_mat& rho, const kraus_ops& ops) {
    int dim = rho.n_rows;
    int n_qubits = slog2(dim);
    cx_mat result = cx_mat(dim, dim, zeros);
    int n_ops = ops.size();

    int total_combinations = pow(n_ops, n_qubits);
    for (int comb = 0; comb < total_combinations; comb++) {
        kraus_ops ops_for_qubits;
        int temp = comb;

        // Determine which Kraus operator to use for each qubit
        for (int q = 0; q < n_qubits; q++) {
            int op_idx = temp % n_ops;
            temp /= n_ops;
            ops_for_qubits.push_back(ops[op_idx]);
        }

        // Build full Kraus operator as tensor product
        if (ops_for_qubits.empty()) return cx_mat(1, 1, ones);

        cx_mat K = ops_for_qubits[0];
        for (size_t i = 1; i < ops_for_qubits.size(); i++) {
            K = kron(K, ops_for_qubits[i]);
        }
        // Apply Kraus operator
        result += (K * rho) * conj(K).t();
    }
    return result;
}
