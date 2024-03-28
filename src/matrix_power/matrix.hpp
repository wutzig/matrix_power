#include <cstdint>
#include <vector>

extern "C" {
int dgemm_(
    char *transa, char *transb,
    long *m, long *n, long *k,
    double *alpha,
    double *a, long *lda, 
    double *b, long *ldb,
    double *beta,
    double *c, long *ldc
);
int dgeev_(
    char *jobvl, char *jobvr,
    long *n,
    double *a, long *lda,
    double *wr, double *wi,
    double *vl, long *ldvl,
    double *vr, long *ldvr,
    double *work, long *lwork,
    long *info
);
}

namespace matrix_power {
class Matrix {
public:
    Matrix() = default;
    Matrix(uint32_t);
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;
    Matrix(Matrix&&);
    Matrix& operator=(Matrix&&);
    ~Matrix();
    
    Matrix& copy_from(const Matrix&);
    Matrix times(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    void scale(double);
    double max_abs() const;
    void print() const;
    Matrix svd(std::vector<double>&);

private:
    int m_dimension;
    double* m_data = nullptr;
};
}