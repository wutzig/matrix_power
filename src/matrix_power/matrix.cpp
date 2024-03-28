
#include <algorithm>
#include <iostream>
#include <matrix_power/matrix.hpp>

namespace matrix_power {
Matrix::Matrix(uint32_t dim) : m_dimension(dim) {
    m_data = new double[m_dimension * m_dimension];
    for(uint32_t i = 0; i < m_dimension * m_dimension; i++) {
        m_data[i] = i + 1;
    }
}
Matrix::Matrix(Matrix&& matrix) : m_dimension(matrix.m_dimension) {
    m_data = matrix.m_data;
    matrix.m_data = nullptr;
}
Matrix& Matrix::copy_from(const Matrix& matrix) {
    m_dimension = matrix.m_dimension;
    delete[] m_data;
    m_data = new double[m_dimension * m_dimension];
    std::copy(
        matrix.m_data,
        matrix.m_data + (matrix.m_dimension * matrix.m_dimension),
        m_data
    );
    return *this;
}
Matrix& Matrix::operator=(Matrix&& matrix) {
    m_dimension = matrix.m_dimension;
    m_data = matrix.m_data;
    matrix.m_data = nullptr;
    return *this;
}
Matrix Matrix::times(const Matrix& matrix) {
    Matrix answer(m_dimension);
    
    char transa ='N';
    char transb ='N';
    long dimension = static_cast<long>(m_dimension);
    double alpha = 1.0;
    double beta = 0;
    
    dgemm_(
        &transb, &transa,
        &dimension, &dimension, &dimension,
        &alpha,
        matrix.m_data, &dimension,
        m_data, &dimension,
        &beta, 
        answer.m_data, &dimension
    );

    return answer;
}
void Matrix::print() {
    for(uint32_t i = 0; i < m_dimension; i++) {
        for(uint32_t j = 0; j < m_dimension; j++) {
            std::cout << m_data[m_dimension * i + j] << '\t';
        }
        std::cout << std::endl;
    }
}
Matrix::~Matrix() {
    delete[] m_data;
}
}