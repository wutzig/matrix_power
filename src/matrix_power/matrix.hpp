#include <cstdint>

namespace matrix_power {
extern "C" {
int dgemm_(
    char *transa, char *transb, long *m, long *n,
    long *k, double *alpha, double *a, long *lda, 
    double *b, long *ldb, double *beta, double *c__, 
    long *ldc
);
}
class Matrix {
public:
    Matrix() = default;
    Matrix(uint32_t);
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;
    Matrix(Matrix&&);
    Matrix& operator=(Matrix&&);
    Matrix& copy_from(const Matrix&);
    Matrix times(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    void scale(double);
    double max_abs() const;
    void print() const;
    ~Matrix();

private:
    int m_dimension;
    double* m_data = nullptr;
};
}