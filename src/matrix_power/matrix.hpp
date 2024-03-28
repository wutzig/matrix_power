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
    Matrix(uint32_t dim = 3);
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;
    Matrix(Matrix&&);
    Matrix& operator=(Matrix&&);
    Matrix& copy_from(const Matrix&);
    Matrix times(const Matrix&);
    void print();
    ~Matrix();

private:
    int m_dimension;
    double* m_data;
};
}