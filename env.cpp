#include "env.hpp"

static auto idx = std::unordered_map<cx_mat, size_t, mat_hash, mat_equal>{};
static auto matrices = std::vector<const cx_mat*>{};

struct mat_hash
{
    size_t operator()(const cx_mat& m) const
    {
        auto hs = std::hash<size_t>{};
        auto hd = std::hash<double>{};
        auto res = hs(m.n_rows);
        res = (res << 1) ^ hs(m.n_cols);
        for (const auto& cx : m) {
            res = (res << 1) ^ hd(cx.real());
            res = (res << 1) ^ hd(cx.imag());
        }
        return res;
    }
};

struct mat_equal
{
    bool operator()(const cx_mat& m1, const cx_mat& m2) const
    {
        return (m1.n_rows == m2.n_rows) && (m1.n_cols == m2.n_cols) && std::equal(m1.begin(), m1.end(), m2.begin());
    }
};

static const cx_mat* load(int id)
{
    if (id < 0 || matrices.size() <= static_cast<size_t>(id))
        return nullptr;
    return matrices[id];
}

static int store(const cx_mat& m)
{
    auto [it, inserted] = idx.emplace(std::move(m), matrices.size());
    if (inserted)
        matrices.emplace_back(&it->first);
    return it->second;
}


auto mat_row_minor(dim_t rows, dim_t cols, const cx_double* matrix) -> cx_mat { return cx_mat(matrix, rows, cols); }

auto mat_row_minor(dim_t rows, dim_t cols, const double* matrix) -> cx_mat
{
    return mat_row_minor(rows, cols, (const cx_double*)matrix);
}

auto mat_row_major(dim_t rows, dim_t cols, const cx_double* matrix) -> cx_mat
{
    auto res = cx_mat(rows, cols);
    for (dim_t row = 0; row < rows; ++row)
        for (dim_t col = 0; col < cols; ++col)
            res(row, col) = matrix[row * cols + col];
    return res;
}

auto mat_row_major(dim_t rows, dim_t cols, const double* matrix) -> cx_mat
{
    return mat_row_major(rows, cols, (const cx_double*)matrix);
}

void from_mat(const cx_mat& matrix, dim_t rows, dim_t cols, double* ret)
{
    assert(matrix.n_rows == rows);
    assert(matrix.n_cols == cols);
    const auto* result = matrix.memptr();
    for (auto i = 0u; i < rows * cols; ++i) {
        ret[2 * i] = result[i].real();
        ret[2 * i + 1] = result[i].imag();
    }
}

void from_mat_row_major(const cx_mat& matrix, dim_t rows, dim_t cols, double* ret)
{
    assert(matrix.n_rows == rows);
    assert(matrix.n_cols == cols);
    for (auto row = 0u; row < rows; ++row)
        for (auto col = 0u; col < cols; ++col) {
            ret[row * cols * col * 2] = matrix(row, col).real();
            ret[row * cols * col * 2 + 1] = matrix(row, col).imag();
        }
}
