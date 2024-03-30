#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <matrix_power/matrix.hpp>

matrix_power::Matrix power_for(matrix_power::Matrix&& matrix, uint32_t power) {
    matrix_power::Matrix orig;
    orig.copy_from(matrix);
    for(int32_t p = 0; p < power - 1; p++) {
        matrix = matrix.times(orig);
    }
    return std::move(matrix);
}

matrix_power::Matrix power_recursive(matrix_power::Matrix&& matrix, const matrix_power::Matrix& orig, uint32_t power) {
    if (power == 1) {
        return std::move(matrix);
    }
    
    uint32_t half_power = std::floor(power / 2);
    matrix = power_recursive(std::move(matrix), orig, half_power);
    matrix = matrix.times(matrix);
    if(power % 2 == 1) {
        matrix = matrix.times(orig);
    }
    return std::move(matrix);
}

matrix_power::Matrix power_wtf(matrix_power::Matrix&& matrix, uint32_t power) {
    matrix_power::Matrix answer = matrix_power::Matrix::make_unit(matrix.get_dimension());
    while(power > 0) {
        if(power % 2 == 1) {
            answer = answer.times(matrix);
        }
        matrix = matrix.times(matrix);
        power = power >> 1;
    }
    return answer;
}

int main()
{
    using Matrix = matrix_power::Matrix;
    
    constexpr size_t n_dim = 3;
    constexpr size_t n_pow = 3;
    
    uint32_t dimensions[n_dim] {10, 100, 1000};
    uint32_t powers[n_pow] {10, 100, 1000};

    for(int dim_i = 0; dim_i < n_dim; dim_i++)
    {
        std::cout << "Dimension " << dimensions[dim_i] << std::endl; 
        for(int pow_i = 0; pow_i < n_pow; pow_i++)
        {
            std::cout << "Power " << powers[pow_i] << std::endl;
            Matrix matrix_for(dimensions[dim_i]);
            
            auto t1 = std::chrono::high_resolution_clock::now();
            matrix_for = power_for(std::move(matrix_for), powers[pow_i]);
            auto t2 = std::chrono::high_resolution_clock::now();
            
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms ";
            
            Matrix matrix_recursive(dimensions[dim_i]);
            t1 = std::chrono::high_resolution_clock::now();
            Matrix original;
            original.copy_from(matrix_recursive);
            matrix_recursive = power_recursive(std::move(matrix_recursive), original, powers[pow_i]);
            t2 = std::chrono::high_resolution_clock::now();
            
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms ";
            
            Matrix matrix_wtf(dimensions[dim_i]);
            t1 = std::chrono::high_resolution_clock::now();
            matrix_wtf = power_wtf(std::move(matrix_wtf), powers[pow_i]);
            t2 = std::chrono::high_resolution_clock::now();
            
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms ";
            
            Matrix matrix_svd(dimensions[dim_i]);
            t1 = std::chrono::high_resolution_clock::now();
            matrix_svd = matrix_svd.power(powers[pow_i]);
            t2 = std::chrono::high_resolution_clock::now();
            
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms\n";
            std::cout << "\tDifference rec " << (matrix_for - matrix_recursive).max_abs() << std::endl;
            std::cout << "\tDifference wtf " << (matrix_for - matrix_wtf).max_abs() << std::endl;
            std::cout << "\tDifference svd " << (matrix_for - matrix_svd).max_abs() << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}