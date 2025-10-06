#include "dmqs.h"
#include "gates.h"
/// @brief Takes a binary string and converts it into a density matrix with the basis state 01
/// @param bin A binary string e.g "0101"
/// @return The density matrix 
cx_mat BinaryStringToDensityMatrix(const string bin){
    int n = bin.length();

    cx_mat density_matrix;
    if (bin[0] == '0') {
        density_matrix = S0;
    } else {
        density_matrix = S1;
    }

    for (int i = 1; i < n; i++) {
        if (bin[i] == '0') {
            density_matrix = kron(density_matrix, S0);
        } else {
            density_matrix = kron(density_matrix, S1);
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
    cx_mat UN = target == 0 ? U1 : Id;
    for(int i = 1; i < n; i++){
        if (i == target){
            UN = kron(UN, U1);
        } else {
            UN = kron(UN, Id);
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

/// @brief Checks equality between two density matrixes
/// @param rho1 First density matrix
/// @param rho2 Second density matrix
/// @param delta Maximum difference allowed for equality 
/// @return Whether the two density matrixes are equal
bool DensityMatrixApproxEq(cx_mat rho1, cx_mat rho2, double delta) {
    assert(rho1.n_rows == rho2.n_rows);
    assert(rho1.n_cols == rho2.n_cols);
    assert(rho1.n_cols == rho2.n_rows);
    for(int i = 0; i < rho1.n_rows; i++) {
        for(int j = 0; j < rho1.n_cols; j++){
            if ((rho1.at(i,j).real() - rho2.at(i,j).real()) > delta ) {
                return false;
            } 
            if ((rho1.at(i,j).imag() - rho2.at(i,j).imag()) > delta ) {
                return false;
            } 
        }
    }
    return true;
}

int rearrangeBits(int i, vector<int> a) {
    int ret = 0;
    for (int j = 0; j < a.size(); j++){
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
    assert(rho.n_rows == rho.n_cols);
    assert(rho.n_rows >= targets.size());
    int n = ceil(log2(rho.n_rows));
    int traced = 1 << (n-targets.size());
    int keept = 1 << targets.size();
    vector<int> tt;
    for(int i = 0; i < n; i++){
        tt.push_back(i);
    }
    int offset = 0;
    for (int i = 0; i < targets.size(); i++) {
        tt.erase(tt.begin() + targets[i] + offset);
        offset -= 1;
    }
    cx_mat result = cx_mat(keept,keept, fill::zeros);
    for (int bit = 0; bit < traced; bit++) {
        int shbr = rearrangeBits(bit, tt);
        for(int r = 0; r < traced; r++) {
            int ir = shbr | rearrangeBits(r, targets);
            for(int c = 0; c < traced; c++){
                int ic = shbr | rearrangeBits(c, targets);
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
vector<int> Sample(cx_mat rho, vector<double> random_values) {
    vector<int> samples;
    vector<double> weights = GetPropabilities(rho);
    samples.reserve(random_values.size());
    
    // Precompute cumulative distribution once
    vector<double> cumulative(weights.size());
    partial_sum(weights.begin(), weights.end(), cumulative.begin());
    double total = cumulative.back();
    for (double& val : cumulative) {
        val /= total;
    }
    
    // Sample for each random value
    for (double r : random_values) {
        auto it = lower_bound(cumulative.begin(), cumulative.end(), r);
        samples.push_back(distance(cumulative.begin(), it));
    }
    
    return samples;
}