#include "dmqs/channels.hpp"
#include <vector>
// Amplitude damping channel Kraus operators
vector<cx_mat> amplitude_damping_ops(double p) {
    vector<cx_mat> ops;
    auto E0 = cx_mat::fixed<2, 2>(zeros);
    E0(0, 0) = 1;
    E0(1, 1) = sqrt(1 - p);
    auto E1 = cx_mat::fixed<2, 2>(zeros);
    E1(0, 1) = sqrt(p);
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

// Phase damping channel Kraus operators
vector<cx_mat> phase_damping_ops(double p) {
    vector<cx_mat> ops;
    auto E0 = cx_mat::fixed<2, 2>(zeros);
    E0(0, 0) = 1;
    E0(1, 1) = sqrt(1 - p);
    auto E1 = cx_mat::fixed<2, 2>(zeros);
    E1(1, 1) = sqrt(p);
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

// Bit flip channel Kraus operators
vector<cx_mat> bit_flip_ops(double p) {
    vector<cx_mat> ops;
    auto E0 = cx_mat::fixed<2, 2>(zeros);
    E0(0, 1) = 1 * sqrt(1 - p);
    E0(1, 0) = 1 * sqrt(1 - p);
    auto E1 = cx_mat::fixed<2, 2>(zeros);
    E1(0, 0) = sqrt(p);
    E1(1, 1) = sqrt(p);
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

vector<cx_mat> phase_flip_ops(double p) {
    vector<cx_mat> ops;
    auto E0 = cx_mat::fixed<2, 2>(zeros);
    E0(0, 0) = sqrt(p);
    E0(1, 1) = sqrt(p);
    auto E1 = cx_mat::fixed<2, 2>(zeros);
    E1(0, 0) = sqrt(1 - p);
    E1(1, 1) = sqrt(1 - p);
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}

// Bit Phase flip channel Kraus operators
vector<cx_mat> bit_phase_flip_ops(double p) {
    vector<cx_mat> ops;
    auto E0 = cx_mat::fixed<2, 2>(zeros);
    E0(0, 0) = sqrt(p);
    E0(1, 1) = sqrt(p);
    auto E1 = cx_mat::fixed<2, 2>(zeros);
    E1(0, 1) = cx_double(0, -sqrt(1 - p));
    E1(1, 0) = cx_double(0, sqrt(1 - p));
    ops.push_back(E0);
    ops.push_back(E1);
    return ops;
}


vector<cx_mat> depolarizing_ops(double p) {
    vector<cx_mat> ops;
    double sqrt_p_over_3 = sqrt(p / 3.0);
    auto E0 = cx_mat::fixed<2, 2>(zeros);
    E0(0, 0) = sqrt(1 - p);
    E0(1, 1) = sqrt(1 - p);
    auto E1 = cx_mat::fixed<2, 2>(zeros);
    E1(0, 1) = sqrt_p_over_3;
    E1(1, 0) = sqrt_p_over_3;

    auto E2 = cx_mat::fixed<2, 2>(zeros);
    E2(0, 1) = cx_double(0, -sqrt_p_over_3);
    E2(1, 0) = cx_double(0, sqrt_p_over_3);
    auto E3 = cx_mat::fixed<2, 2>(zeros);
    E3(0, 0) = sqrt_p_over_3;
    E3(1, 1) = -sqrt_p_over_3;
    ops.push_back(E0);
    ops.push_back(E1);
    ops.push_back(E2);
    ops.push_back(E3);
    return ops;
}

std::function<vector<cx_mat>(double)> u_channel_to_ops_f(u_channel channel) {
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

cx_mat apply_channel(const cx_mat& rho, const vector<cx_mat>& kraus_ops) {
    int dim = rho.n_rows;
    int n_qubits = slog2(dim);
    cx_mat result = cx_mat(dim, dim, zeros);
    int n_ops = kraus_ops.size();

    int total_combinations = pow(n_ops, n_qubits);

    for (int comb = 0; comb < total_combinations; comb++) {
        vector<cx_mat> ops_for_qubits;
        int temp = comb;

        // Determine which Kraus operator to use for each qubit
        for (int q = 0; q < n_qubits; q++) {
            int op_idx = temp % n_ops;
            temp /= n_ops;
            ops_for_qubits.push_back(kraus_ops[op_idx]);
        }

        // Build full Kraus operator as tensor product
        if (ops_for_qubits.empty()) return cx_mat(1, 1, ones);

        cx_mat sum = cx_mat(2, 2, zeros);
        for (auto op : ops_for_qubits) {
            sum += op * adjoint(op);
        }

        cx_mat K = ops_for_qubits[0];
        for (size_t i = 1; i < ops_for_qubits.size(); i++) {
            K = kron(K, ops_for_qubits[i]);
        }
        // Apply Kraus operator
        result += (K * rho) * conj(K).t();
    }
    return result;
}

cx_mat apply_channel_with_ops_for_each_qubit(const cx_mat& rho,
                                             const vector<cx_mat>& kraus_ops) {
    int dim = rho.n_rows;
    int n_qubits = slog2(dim);
    cx_mat result = cx_mat(dim, dim, zeros);
    int n_ops_qb = kraus_ops.size() / n_qubits;

    int total_combinations = pow(n_ops_qb, n_qubits);
    for (int comb = 0; comb < total_combinations; comb++) {
        vector<cx_mat> ops_for_qubits;
        int temp = comb;

        // Determine which Kraus operator to use for each qubit
        for (int q = 0; q < n_qubits; q++) {
            int op_idx = temp % (n_ops_qb);
            temp /= (n_ops_qb);
            ops_for_qubits.push_back(kraus_ops[q * (n_ops_qb) + op_idx]);
        }

        // Build full Kraus operator as tensor product
        if (ops_for_qubits.empty()) return cx_mat(1, 1, ones);

        cx_mat K_full = ops_for_qubits[0];
        for (size_t i = 1; i < ops_for_qubits.size(); i++) {
            K_full = kron(K_full, ops_for_qubits[i]);
        }

        // Apply Kraus operator
        result += K_full * rho * conj(K_full).t();
    }
    return result;
}
