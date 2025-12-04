#include <dmqs/channels.hpp>
#include <vector>
// Amplitude damping channel Kraus operators
vector<kraus_t> amplitude_damping_ops(const double& p) {
    return {
        { {cx_double(1, 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(1 - p), 0)} },
        { {cx_double(0, 0), cx_double(sqrt(p), 0)},
          {cx_double(0, 0), cx_double(0, 0)}}
    };
}

vector<kraus_t> amplitude_damping_ops2(const double& p, const double& t) {
  return {
        { {cx_double(1, 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(exp(-p*t)), 0)} },
        { {cx_double(0, 0), cx_double(sqrt(1 - exp(-p*t)), 0)},
          {cx_double(0, 0), cx_double(0, 0)}}
    };
}

vector<kraus_t> generalized_amplitude_damping_ops(const double& p,
                                            const double& gamma) {
    return {
        { { cx_double(sqrt(p), 0), cx_double(0, 0) },
          { cx_double(0, 0), cx_double(sqrt(p) * sqrt(1 - gamma), 0) } },
        { { cx_double(0, 0), cx_double(sqrt(p) * sqrt(gamma), 0) },
          { cx_double(0, 0), cx_double(0, 0) } },
        { { cx_double(sqrt(1 - p) * sqrt(1 - gamma), 0), cx_double(0, 0) },
          { cx_double(0, 0), cx_double(sqrt(1 - p), 0) } },
        { { cx_double(0, 0), cx_double(0, 0) },
          {cx_double(sqrt(1 - p) * sqrt(gamma), 0), cx_double(0, 0) } }
    };
}

// Phase damping channel Kraus operators
vector<kraus_t> phase_damping_ops(const double& p) {
    return {
        { {cx_double(1, 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(1-p), 0)} },
        { {cx_double(0, 0), cx_double(0, 0)},
          {cx_double(0, 0), cx_double(sqrt(p), 0)} },
    };
}

// Bit flip channel Kraus operators
vector<kraus_t> bit_flip_ops(const double& p) {
    return {
        { { cx_double(sqrt(p), 0), cx_double(0, 0) },
          { cx_double(0, 0), cx_double(sqrt(p), 0) } },
        { { cx_double(0, 0), cx_double(sqrt(1 - p), 0) },
          { cx_double(sqrt(1 - p), 0), cx_double(0, 0) } }
    };
}

vector<kraus_t> phase_flip_ops(const double& p) {
    return {
        { { cx_double(sqrt(p), 0), cx_double(0, 0) },
          { cx_double(0, 0), cx_double(sqrt(p), 0) } },
        { { cx_double(sqrt(1 - p), 0), cx_double(0, 0) },
          { cx_double(0, 0), cx_double(sqrt(1 - p)) } }
    };
}

// Bit Phase flip channel Kraus operators
vector<kraus_t> bit_phase_flip_ops(const double& p) {
    return {
        { { cx_double(sqrt(p), 0), cx_double(0, 0) },
          { cx_double(0, 0), cx_double(sqrt(p), 0) } },
        { { cx_double(0, 0), cx_double(0, -sqrt(1 - p)) },
          { cx_double(0, sqrt(1 - p)), cx_double(0, 0) } }
    };
}

vector<kraus_t> depolarizing_ops(const double& p) {
    double sqrt1 = sqrt(1 - ((p * 3.0)/ 4.0));
    double sqrt2 = sqrt(p / 4.0);
    return {
      { { cx_double(sqrt1, 0), cx_double(0, 0) },
        { cx_double(0, 0), cx_double(sqrt1, 0) } },
      { { cx_double(0, 0), cx_double(sqrt2, 0) },
        { cx_double(sqrt2, 0), cx_double(0, 0) } },
      { { cx_double(0, 0), cx_double(0, -sqrt2) },
        { cx_double(0, sqrt2), cx_double(0, 0) } },
      { { cx_double(sqrt2, 0), cx_double(0, 0) },
        { cx_double(0, 0), cx_double(-sqrt2, 0) } },
    };
}

channel_t u_channel_to_ops_f(u_channel channel) {
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

cx_mat apply_channel(const cx_mat& rho, const vector<kraus_t>& ops) {
    int dim = rho.n_rows;
    int n_qubits = slog2(dim);
    cx_mat result = cx_mat(dim, dim, zeros);
    int n_ops = ops.size();

    int total_combinations = pow(n_ops, n_qubits);
    for (int comb = 0; comb < total_combinations; comb++) {
        vector<kraus_t> ops_for_qubits;
        int temp = comb;
        //std::cout << "Run: " << comb << std::endl;
        // Determine which Kraus operator to use for each qubit
        for (int q = 0; q < n_qubits; q++) {
            int op_idx = temp % n_ops;
            //std::cout << op_idx << "," << std::endl;
            temp /= n_ops;
            ops_for_qubits.push_back(ops[op_idx]);
        }

        cx_mat K = ops_for_qubits[0];
        for (size_t i = 1; i < ops_for_qubits.size(); i++) {
            K = kron(K, ops_for_qubits[i]);
        }
        // Apply Kraus operator
        result += (K * rho) * K.t();
    }
    return result;
}
