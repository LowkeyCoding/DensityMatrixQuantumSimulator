#include "dmqs.hpp"

/// @brief Takes a binary string and converts it size_to a density matrix with the basis state 01
/// @param bin A binary string e.g "0101"
/// @return The density matrix 
cx_mat BinaryStringToDensityMatrix(const string bin){
    size_t n = bin.length();

    cx_mat density_matrix;
    if (bin[0] == '0') {
        density_matrix = B0();
    } else {
        density_matrix = B1();
    }

    for (size_t i = 1; i < n; i++) {
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
cx_mat GateToNQubitSystem(cx_mat U1, size_t target, size_t n){
    cx_mat UN = target == 0 ? U1 : Id();
    for(size_t i = 1; i < n; i++){
        if (i == target){
            UN = kron(UN, U1);
        } else {
            UN = kron(UN, Id());
        }
    }
    return UN;
}

/// @brief Applies a Gate U to the density matrix rho
/// @param rho 
/// @param U 
/// @return The modified density matrix
cx_mat ApplyGateToDensityMatrix(cx_mat rho, cx_mat U){
    assert(rho.n_rows == U.n_rows);
    assert(rho.n_cols == U.n_cols);
    assert(rho.n_cols == U.n_rows);
    return (U * rho) * (conj(U).t());
}

bool IsPure(cx_mat rho, double delta) {
    cx_double t = trace(rho*rho);
    return t.real() - 1 < delta;
}


cx_mat ApplyGate(cx_mat rho, u_gate gate, size_t qubit) {
    size_t qubit_count = ceil(log2(rho.n_rows));
    cx_mat U;
    switch (gate) {
        case GX:
            U = X();
            break;
        case GY:
            U = Y();
            break;
        case GZ:
            U = Z();
            break;
        case GCX:
            return ApplyGateToDensityMatrix(rho, CX());
            break;
        case GH:
            U = H();
            break;
        case GB0:
            U = B0();
            break;
        case GB1:
            U = B1();
            break;
        default:
            U = Id();
            break;
    }
    return ApplyGateToDensityMatrix(rho, GateToNQubitSystem(U,qubit,qubit_count));
}

size_t rearrangeBits(size_t i, vector<size_t> a) {
    size_t ret = 0;
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
cx_mat PartialTrace(cx_mat rho, vector<size_t> targets) {
    assert(rho.n_rows == rho.n_cols);
    assert(rho.n_rows >= targets.size());
    size_t n = ceil(log2(rho.n_rows));
    size_t traced = 1 << (n-targets.size());
    size_t keept = 1 << targets.size();
    vector<size_t> tt;
    for(size_t i = 0; i < n; i++){
        tt.push_back(i);
    }
    size_t offset = 0;
    for (size_t i = 0; i < targets.size(); i++) {
        tt.erase(tt.begin() + targets[i] + offset);
        offset -= 1;
    }
    cx_mat result = cx_mat(keept,keept, fill::zeros);
    for (size_t bit = 0; bit < traced; bit++) {
        size_t shbr = rearrangeBits(bit, tt);
        for(size_t r = 0; r < traced; r++) {
            size_t ir = shbr | rearrangeBits(r, targets);
            for(size_t c = 0; c < traced; c++){
                size_t ic = shbr | rearrangeBits(c, targets);
                result(r,c) += rho(ir,ic);
            }
        }
    }
    return result;
}

/// @brief Gets the weights from the diagonal of the density matrix
/// @param rho Density matrix to get weights from
/// @return The weights of density matrix.
vector<double> GetPropabilities(cx_mat rho) {
    // Get the diagonal elements (these are the probabilities)
    cx_vec diagonal = rho.diag();
    
    // Convert to real probabilities (diagonal of density matrix should be real and non-negative)
    vector<double> weights(diagonal.n_elem);
    for (size_t i = 0; i < diagonal.n_elem; ++i) {
        weights[i] = diagonal(i).real();  // Take real part (imaginary should be zero)
    }
    
    return weights;
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
        throw invalid_argument("Random value must be in [0, 1)");
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