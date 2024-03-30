
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <matrix_power/matrix.hpp>

namespace matrix_power {
Matrix::Matrix(uint32_t dim) : m_dimension(dim) {
    m_data = new double[m_dimension * m_dimension];
    std::fill(m_data, m_data + (m_dimension * m_dimension), 0.0);
    double diag = 10.0 / m_dimension;
    double off_diag = -1.0 / m_dimension;
    m_data[0] = diag;
    m_data[1] = off_diag;
    for(uint32_t i = 1; i < m_dimension - 1; i++) {
        m_data[i * m_dimension + i - 1] = off_diag;
        m_data[i * m_dimension + i] = diag;
        m_data[i * m_dimension + i + 1] = off_diag;
    }
    m_data[m_dimension * m_dimension - 2] = off_diag;
    m_data[m_dimension * m_dimension - 1] = diag;
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
    delete[] m_data;
    m_dimension = matrix.m_dimension;
    m_data = matrix.m_data;
    matrix.m_data = nullptr;
    return *this;
}
Matrix Matrix::times(const Matrix& matrix, bool transpose_a) const {
    Matrix answer(m_dimension);
    
    char transa = transpose_a ? 'T' : 'N';
    char transb = 'N';
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
void Matrix::print() const {
    for(uint32_t i = 0; i < m_dimension; i++) {
        for(uint32_t j = 0; j < m_dimension; j++) {
            std::cout << m_data[m_dimension * i + j] << '\t';
        }
        std::cout << std::endl;
    }
}

Matrix Matrix::operator-(const Matrix& matrix) const {
    Matrix answer(m_dimension);
    for(uint32_t i = 0; i < m_dimension * m_dimension; i++) {
        answer.m_data[i] = m_data[i] - matrix.m_data[i];
    }
    return answer;
}

double Matrix::max_abs() const {
    double answer = std::abs(m_data[0]);
    for(uint32_t i = 0; i < m_dimension * m_dimension; i++) {
        answer = std::max(answer, std::abs(m_data[i]));
    }
    return answer;
}

void Matrix::scale(double factor) {
    for(uint32_t i = 0; i < m_dimension * m_dimension; i++) {
        m_data[i] *= factor;
    }
}

Matrix Matrix::svd(std::vector<double>& eigen_values) const {
    Matrix answer(m_dimension);
    answer.scale(0);

    char jobvl = 'N';
    char jobvr = 'V';
    long dimension = static_cast<long>(m_dimension);
    long lwork = 4 * dimension;
    long info = 0;
    std::vector<double> dummy(m_dimension);
    std::vector<double> work(lwork);

    dgeev_(
        &jobvl, &jobvr,
        &dimension,
        m_data, &dimension,
        eigen_values.data(), dummy.data(),
        nullptr, &dimension,
        answer.m_data, &dimension,
        work.data(), &lwork,
        &info
    );
    return answer;
}

Matrix Matrix::power(uint32_t power) const {
    std::vector<double> eigen_values(m_dimension);

    Matrix eigen_vectors = svd(eigen_values);
    Matrix answer(m_dimension);
    answer.scale(0);
    for(uint32_t i = 0; i < m_dimension; i++) {
        answer.m_data[i * m_dimension + i] = std::pow(eigen_values[i], power);
    }
    answer = eigen_vectors.times(answer, true);
    answer = answer.times(eigen_vectors);

    return answer;
}

Matrix Matrix::make_unit(uint32_t dimension) {
    Matrix answer;
    answer.m_dimension = dimension;
    answer.m_data = new double[answer.m_dimension * answer.m_dimension];
    std::fill(answer.m_data, answer.m_data + (answer.m_dimension * answer.m_dimension), 0.0);
    for(uint32_t i = 0; i < answer.m_dimension; i++) {
        answer.m_data[i * answer.m_dimension + i] = 1.0;
    }
    return answer;
}

uint32_t Matrix::get_dimension() const {
    return m_dimension;
}
Matrix::~Matrix() {
    delete[] m_data;
}
}