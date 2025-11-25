#include <uppaal/uppaal.h>
#include <vector>

// Initialize density matrix with binary state string
// (e.g., "01" for |01⟩ or "+-" for |+-⟩)
// rho_size = 1 << 2*N+1, bin = binary state string of length N
extern "C" void UInitBinState(double* rho, int rho_size, const char* state) {
    assert(size_t(rho_size) == strlen(state));
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat res = BinaryStringToDensityMatrix(state);
    memmove(rho, res.memptr(), mat_size * sizeof(double));
}

// Apply single-qubit gate to target qubit
extern "C" void UApplyGate(double* rho, int rho_size, int gate, int target) {
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    cx_mat res = ApplyGate(in_mat, static_cast<u_gate>(gate), target);
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Apply controlled gate (control -> target)
extern "C" void UApplyCGate(double* rho, int rho_size, int gate, int control,
                            int target) {
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    cx_mat res = ApplyCGate(in_mat, static_cast<u_gate>(gate), control, target);
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Apply custom unitary matrix U (u_size must equal rho_size)
extern "C" void UApplyUnitary(double* rho, int rho_size, double* U,
                              int u_size) {
    if (rho_size != u_size) {
        throw invalid_argument(
            "Density matrix size should be equivilant to unitary. Got " +
            to_string(rho_size) + " and " + to_string(u_size));
    }
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_rho = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    cx_mat in_U = cx_mat(reinterpret_cast<cx_double*>(U),
                         mat_row, mat_row, false, true);
    cx_mat res = ApplyGateToDensityMatrix(in_rho, in_U);
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Partial trace: extract subsystem state into prho
// prho_size = 1 << N*2-1 where N = len(targets)
extern "C" void UPartialTrace(double* rho, int rho_size, double* prho,
                              int prho_size, int* targets, int targets_size) {
    if (prho_size != targets_size) {
        throw invalid_argument(
            "The partial density matrix should have the same "
            "size as the number of targets. Got " +
            to_string(rho_size) + " and " + to_string(targets_size));
    }
    size_t mat_row = 1 << rho_size;
    vector<int> t = vector<int>(targets, targets + targets_size);
    cx_mat in_rho = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    cx_mat res = PartialTrace(in_rho, t);
    size_t mat_size = ((1 << prho_size) * (1 << prho_size)) << 1;
    memcpy(prho, res.memptr(), mat_size * sizeof(double));
}

// Project target qubit onto basis state (|0⟩ or |1⟩)
extern "C" void UBasisProjection(double* rho, int rho_size, int target,
                                 int state) {
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, true, true);
    cx_mat res = BasisProjection(in_mat, target, state);
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Project multiple target qubits onto basis state (bitmask)
extern "C" void UBasisProjections(double* rho, int rho_size, int* targets,
                                  int targets_size, int state) {
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    vector<int> t = vector<int>(targets, targets + targets_size);
    cx_mat res = BasisProjections(in_mat, t, state);
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Measure target qubits and return result (does not collapse state)
// Use UBasisProjections to collapse after measurement
extern "C" int UPartialMeasure(double* rho, int rho_size, int* targets,
                               int targets_size, double r) {
    size_t mat_row = 1 << rho_size;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    vector<int> t = vector<int>(targets, targets + targets_size);
    return PartialSample(in_mat, t, r);
}

// Measure all qubits and return result (does not collapse state)
extern "C" int UMeasureAll(double* rho, int rho_size, double r) {
    size_t mat_row = 1 << rho_size;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    return Sample(in_mat, r);
}

// Apply amplitude damping and dephasing noise channel on rho as seen
// in: 10.1098/rspa.2008.0439
extern "C" void UAmplitudeDampeningAndDephasing(double* rho, int rho_size,
                                                const double* T1,
                                                const double* T2, double t) {
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    cx_mat res = ApplyAmplitudeDampeningAndDephasing(in_mat, T1, T2, t);
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Apply a noise channel to the density matrix.
// See https://link.springer.com/content/pdf/10.1007/s10773-019-04332-z.pdf
// and https://docs.pennylane.ai/en/stable/_modules/pennylane/ops/channel.html
// for details on the implementation of channels.
extern "C" void UApplyChannel(double* rho, int rho_size, int channel,
                              double* probs, int probs_size) {
    auto chan = static_cast<u_channel>(channel);
    vector<cx_mat> ops;
    cx_mat res;
    std::function<vector<cx_mat>(double)> chan_f = nullptr;
    size_t mat_row = 1 << rho_size;
    size_t mat_size = mat_row * mat_row * 2;
    cx_mat in_mat = cx_mat(reinterpret_cast<cx_double*>(rho),
                           mat_row, mat_row, false, true);
    switch (chan) {
        case AMPLITUDE_DAMPING:
            chan_f = amplitude_damping_ops;
            break;
        case PHASE_DAMPING:
            chan_f = phase_damping_ops;
            break;
        case DEPOLARIZING:
            chan_f = depolarizing_ops;
            break;
        case PHASE_FLIP:
            chan_f = phase_flip_ops;
            break;
        case BIT_PHASE_FLIP:
            chan_f = bit_phase_flip_ops;
            break;
        default:
            throw invalid_argument("Unknown channel type");
    }
    if (probs_size > 1) {
        for (int i = 0; i < probs_size; i++) {
            auto ops_for_qubit = chan_f(probs[i]);
            ops.insert(ops.end(), ops_for_qubit.begin(), ops_for_qubit.end());
        }
        res = apply_channel(in_mat, ops);
    } else {
        ops = chan_f(probs[0]);
        res = apply_channel(in_mat, ops);
    }
    memcpy(rho, res.memptr(), mat_size * sizeof(double));
}

// Applies noise with a single probability parameter
extern "C" void UApplySChannel(double* rho, int rho_size, int channel,
                              double probs) {
    UApplyChannel(rho, rho_size, channel, &probs, 1);
}
