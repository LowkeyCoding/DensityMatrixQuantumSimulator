#include "dmqs/channels.hpp"

// Amplitude damping channel Kraus operators
vector<cx_mat> amplitude_damping_ops(double gamma) {
    vector<cx_mat> ops;
    
    cx_mat E0(2, 2, fill::zeros);
    E0(0, 0) = 1.0;
    E0(1, 1) = sqrt(1 - gamma);
    
    cx_mat E1(2, 2, fill::zeros);
    E1(0, 1) = sqrt(gamma);
    
    ops.push_back(E0);
    ops.push_back(E1);
    
    return ops;
}

// Phase damping channel Kraus operators  
vector<cx_mat> phase_damping_ops(double gamma) {
    vector<cx_mat> ops;
    
    cx_mat E0(2, 2, fill::zeros);
    E0(0, 0) = 1.0;
    E0(1, 1) = sqrt(1 - gamma);
    
    cx_mat E1(2, 2, fill::zeros);
    E1(1, 1) = sqrt(gamma);
    
    ops.push_back(E0);
    ops.push_back(E1);
    
    return ops;
}

// Depolarizing channel Kraus operators
vector<cx_mat> depolarizing_ops(double p) {
    vector<cx_mat> ops;
    
    cx_mat I(2, 2, fill::zeros);
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    
    cx_mat X(2, 2, fill::zeros);
    X(0, 1) = 1.0;
    X(1, 0) = 1.0;
    
    cx_mat Y(2, 2, fill::zeros);
    Y(0, 1) = -cx_double(0, 1);  // -i
    Y(1, 0) = cx_double(0, 1);   // i
    
    cx_mat Z(2, 2, fill::zeros);
    Z(0, 0) = 1.0;
    Z(1, 1) = -1.0;
    
    ops.push_back(sqrt(1 - 3*p/4) * I);
    ops.push_back(sqrt(p/4) * X);
    ops.push_back(sqrt(p/4) * Y);
    ops.push_back(sqrt(p/4) * Z);
    
    return ops;
}

// Bit flip channel Kraus operators
vector<cx_mat> bit_flip_ops(double p) {
    vector<cx_mat> ops;
    
    cx_mat E0(2, 2, fill::zeros);
    E0(0, 0) = sqrt(1 - p);
    E0(1, 1) = sqrt(1 - p);
    
    cx_mat E1(2, 2, fill::zeros);
    E1(0, 1) = sqrt(p);
    E1(1, 0) = sqrt(p);
    
    ops.push_back(E0);
    ops.push_back(E1);
    
    return ops;
}

// Phase flip channel Kraus operators
vector<cx_mat> phase_flip_ops(double p) {
    vector<cx_mat> ops;
    
    cx_mat E0(2, 2, fill::zeros);
    E0(0, 0) = sqrt(1 - p);
    E0(1, 1) = sqrt(1 - p);
    
    cx_mat E1(2, 2, fill::zeros);
    E1(0, 0) = sqrt(p);
    E1(1, 1) = -sqrt(p);
    
    ops.push_back(E0);
    ops.push_back(E1);
    
    return ops;
}

// Bit-phase flip channel Kraus operators
vector<cx_mat> bit_phase_flip_ops(double p) {
    vector<cx_mat> ops;
    
    cx_mat E0(2, 2, fill::zeros);
    E0(0, 0) = sqrt(1 - p);
    E0(1, 1) = sqrt(1 - p);
    
    cx_mat E1(2, 2, fill::zeros);
    E1(0, 1) = -cx_double(0, 1) * sqrt(p);  // -i√p
    E1(1, 0) = cx_double(0, 1) * sqrt(p);   // i√p
    
    ops.push_back(E0);
    ops.push_back(E1);
    
    return ops;
}

cx_mat apply_channel(const cx_mat& rho, const vector<cx_mat>& kraus_ops, int n_qubits) {
    int dim = 1 << n_qubits;
    cx_mat result = cx_mat(dim, dim, fill::zeros);
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
        if (ops_for_qubits.empty()) return cx_mat(1, 1, fill::ones);
    
        cx_mat K_full = ops_for_qubits[0];
        for (size_t i = 1; i < ops_for_qubits.size(); i++) {
            K_full = kron(K_full, ops_for_qubits[i]);
        }
        
        // Apply Kraus operator
        result += K_full * rho * K_full.t();
    }
    
    return result;
}


// Kronecker product of multiple matrices
cx_mat multi_kron(const vector<cx_mat>& matrices) {
    if (matrices.empty()) return cx_mat(1, 1, fill::ones);
    
    cx_mat result = matrices[0];
    for (size_t i = 1; i < matrices.size(); i++) {
        result = kron(result, matrices[i]);
    }
    return result;
}