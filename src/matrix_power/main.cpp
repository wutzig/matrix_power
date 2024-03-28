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

matrix_power::Matrix power_wtf(matrix_power::Matrix&& matrix, const matrix_power::Matrix& orig, uint32_t power) {
    if (power == 1) {
        return std::move(matrix);
    }
    
    uint32_t half_power = std::floor(power / 2);
    matrix = power_wtf(std::move(matrix), orig, half_power);
    matrix = matrix.times(matrix);
    if(power % 2 == 1) {
        matrix = matrix.times(orig);
    }
    return std::move(matrix);
}

int main()
{
    using Matrix = matrix_power::Matrix;
    
    constexpr size_t n_dim = 3;
    constexpr size_t n_pow = 3;
    
    uint32_t dimensions[n_dim] {10, 100, 1000};
    uint32_t powers[n_pow] {10, 100, 1000};

    for(int dim_i = 0; dim_i < 3; dim_i++)
    {
        std::cout << "Dimension " << dimensions[dim_i] << std::endl; 
        for(int pow_i = 0; pow_i < 3; pow_i++)
        {
            std::cout << "Power " << powers[pow_i] << std::endl;
            Matrix matrix_for(dimensions[dim_i]);
            
            auto t1 = std::chrono::high_resolution_clock::now();
            matrix_for = power_for(std::move(matrix_for), powers[pow_i]);
            auto t2 = std::chrono::high_resolution_clock::now();
            
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms ";
            
            Matrix matrix_wtf(dimensions[dim_i]);
            Matrix original;
            original.copy_from(matrix_wtf);
            t1 = std::chrono::high_resolution_clock::now();
            matrix_wtf = power_wtf(std::move(matrix_wtf), original, powers[pow_i]);
            t2 = std::chrono::high_resolution_clock::now();
            
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms\n";
            std::cout << "\tDifference " << (matrix_for - matrix_wtf).max_abs() << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}