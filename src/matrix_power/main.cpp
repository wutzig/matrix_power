#include <matrix_power/matrix.hpp>
int main()
{
    using Matrix = matrix_power::Matrix;
    Matrix mat(3);
    mat.print();
    mat = mat.times(mat);
    mat.print();
    return 0;
}