#include "dmqs/dmqs.hpp"

/// @brief Takes a binary string and converts it into a density matrix with the basis state 01
/// @param bin A binary string e.g "0101"
/// @return The density matrix 
cx_mat BinaryStringToDensityMatrix(const string bin){
    int n = bin.length();

    cx_mat density_matrix;
    if (bin[0] == '0') {
        density_matrix = B0();
    } else {
        density_matrix = B1();
    }

    for (int i = 1; i < n; i++) {
        if (bin[i] == '0') {
            density_matrix = kron(density_matrix, B0());
        } else {
            density_matrix = kron(density_matrix, B1());
        }
    }
    
    return density_matrix;
}

/// @brief Creates a gate that applies a single-qubit operation to a specific target qubit in an n-qubit system.
/// @param U1 The single-qubit gate to be applied.
/// @param target The index (zero-based) of the target qubit within the system.
/// @param n The total number of qubits in the system.
/// @return An n-qubit gate that only affects the specified target qubit.
cx_mat GateToNQubitSystem(cx_mat U1, int target, int n){
    if (target > 0){
        U1 = kron(Id(target), U1);
    }
    int diff = int(abs(target - n)) - 1;
    if (diff > 0) {
        U1 = kron(U1,Id(diff));
    }
    return U1;
}

/// @brief Applies a Gate U to the density matrix rho
/// @param rho 
/// @param U 
/// @return The modified density matrix
cx_mat ApplyGateToDensityMatrix(cx_mat rho, cx_mat U){
    return (U * rho) * (conj(U).t());
}

bool IsPure(cx_mat rho, double delta) {
    cx_double t = trace(rho*rho);
    return t.real() - 1 < delta;
}


cx_mat ApplyGate(cx_mat rho, u_gate gate, int qubit) {
    int qubit_count = ceil(log2(rho.n_rows));
    cx_mat U;
    switch (gate) {
        case GID:
            U = Id();
            break;
        case GX:
            U = X();
            break;
        case GY:
            U = Y();
            break;
        case GZ:
            U = Z();
            break;
        case GH:
            U = H();
            break;
        default:
            throw invalid_argument("ApplyGate " + to_string((int)gate) + " is an invalid gate");
            break;
    }
    return ApplyGateToDensityMatrix(rho, GateToNQubitSystem(U,qubit,qubit_count));
}

cx_mat ApplyCGate(cx_mat rho, u_gate gate, int target, int control) {
    int qubit_count = ceil(log2(rho.n_rows));
    cx_mat U;
    switch (gate) {
        case GCX:
            U = CG(X(),target,control);
            break;
        case GCY:
            U = CG(Y(),target,control);
            break;
        case GCZ:
            U = CG(Z(),target,control);;
            break;
        case GCH:
            U = CG(H(),target,control);
            break;
        default:
            throw invalid_argument("ApplyCGate " + to_string((int)gate) + " is an invalid gate");
            break;
    }
    if (target > 0){
        U = kron(Id(target), U);
    }
    int diff = int(abs(control - qubit_count)) - 1;
    if (diff > 0) {
        U = kron(U,Id(diff));
    }

    return ApplyGateToDensityMatrix(rho, U);
}

int rearrangeBits(int i, vector<int> a) {
    int ret = 0;
    for (size_t j = 0; j < a.size(); j++){
        if (a[j] >= 0) {
            ret |= ((i >> j) & 1) << a[j];
        }
    }
    return ret;
}

/// @brief Takes a partial trace of the provided target qubits.
/// @param rho Density matrix to take a partial trace from.
/// @param targets A vector containing the qubits to trace.
/// @return A denisty matrix of the traced qubits.
cx_mat PartialTrace(cx_mat rho, vector<int> targets) {
    if (rho.n_rows != rho.n_cols) {
        throw invalid_argument("Matrix provided should be square not " + to_string(rho.n_rows) + " by " + to_string(rho.n_cols));
    }
    auto minmax = minmax_element(targets.begin(), targets.end());
    if (*minmax.first < 0 || (size_t)*minmax.second > rho.n_rows) {
        throw invalid_argument("Targets should be greater than 0 and less than " + to_string(rho.n_rows) + ". Min max " + to_string(*minmax.first) + "," + to_string(*minmax.second));
    }
    int n = ceil(log2(rho.n_rows));
    int traced_size = 1 << (n - targets.size()); // size of kept system
    int kept_size = 1 << targets.size(); // size of traced-out system
    // Convert qubit indices to bit positions (qubit 0 is MSB)
    vector<int> bit_positions;
    for (int i = 0; i < n; i++) {
        bit_positions.push_back(n - 1 - i);
    }
    
    vector<int> kept_qubits;
    for (int i = 0; i < n; i++) {
        kept_qubits.push_back(i);
    }
    // Remove target qubits from kept_qubits
    sort(targets.begin(), targets.end(), greater<int>());
    for (int t : targets) {
        kept_qubits.erase(kept_qubits.begin() + t);
    }
    
    // Convert kept_qubits and targets to bit positions
    vector<int> kept_bits;
    for (int k : kept_qubits) {
        kept_bits.push_back(bit_positions[k]);
    }
    vector<int> target_bits;
    for (int t : targets) {
        target_bits.push_back(bit_positions[t]);
    }
    
    cx_mat result = cx_mat(kept_size, kept_size, fill::zeros);
    
    for (int r = 0; r < kept_size; r++) {
        int r_bits = rearrangeBits(r, target_bits);
        for (int c = 0; c < kept_size; c++) {
            int c_bits = rearrangeBits(c, target_bits);
            for (int bit = 0; bit < traced_size; bit++) {
                int bit_bits = rearrangeBits(bit, kept_bits);
                int ir = bit_bits | r_bits;
                int ic = bit_bits | c_bits;
                result(r, c) += rho(ir, ic);
            }
        }
    }
    return result;
}


cx_mat MeasurementGate(const cx_mat rho, int target, double random_value) {
    // Calculate measurement probabilities directly
    int sample = Sample(rho, random_value);
    auto U = sample ? B1() : B0();
    cx_mat rho_projected = ApplyGateToDensityMatrix(rho, GateToNQubitSystem(U,target,ceil(log2(rho.n_rows))));
    
    // Normalize by trace
    cx_double trace_val = trace(rho_projected);
    rho_projected = rho_projected / trace_val;
    
    return rho_projected;
}

/// @brief Gets a sample from the density matrix for each random value provided.
/// @param rho Density matrix to sample from.
/// @return A vector of samples from the density matrix.
int Sample(const cx_mat rho, double random_value) {
    // Validate input
    if (rho.n_rows != rho.n_cols) {
        throw invalid_argument("Density matrix must be square");
    }
    
    if (random_value < 0.0 || random_value >= 1.0) {
        throw invalid_argument("Random value must be in [0, 1]. " + to_string(random_value) +  "was supplied");
    }
    
    // Diagonalize the density matrix to get eigenvalues (probabilities)
    vec eigval = abs(rho.diag());

    // Ensure eigenvalues are non-negative (probabilities)
    if (any(eigval < -1e-10)) {
        throw invalid_argument("Density matrix has negative eigenvalues");
    }

    // Cumulative distribution for sampling
    vec cumulative = cumsum(eigval);

    // Sample based on the random value
    for (uword i = 0; i < cumulative.n_elem; ++i) {
        if (random_value < cumulative(i)) {
            return static_cast<int>(i);
        }
    }
    
    // Fallback: return last index (shouldn't happen if random_value < 1.0)
    return static_cast<int>(cumulative.n_elem - 1);
}

/// @brief Applies amplitude dampening and dephasing channel as seen in: 10.1098/rspa.2008.0439
/// @param rho Density matrix to apply 
/// @param T1 Energy relaxation time
/// @param T2 Phase choherence time
/// @param t time channel acts upon qubits
/// @return 
cx_mat ApplyAmplitudeDampeningAndDephasing(cx_mat rho, double* T1, double* T2, double t) {
    cx_mat temp_state = rho;
    double n = ceil(log2(rho.n_rows));
    for(int i = 0; i < n; i++){
        double px = (1 - exp(-t/T1[i]))*0.25;
        double py = px;
        double pz = 0.5 - px - (exp(-t/T2[i])*0.5);
        cx_mat XN = GateToNQubitSystem(X(),i,n);
        cx_mat YN = GateToNQubitSystem(Y(),i,n);
        cx_mat ZN = GateToNQubitSystem(Z(),i,n);
        temp_state = ((1 - px - py - pz) * temp_state) + (px * XN * temp_state * XN)  + (py * YN * temp_state * YN) + (pz * ZN * temp_state * ZN);
    }
    return temp_state;
}