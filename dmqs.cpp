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
bool DensitryMatrixApproxEq(cx_mat rho1, cx_mat rho2, double delta) {
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