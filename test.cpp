#include "HEADERS.h"
#include "program.cpp"
using namespace LinearAlgebra;

int main(int argc, char *argv[])
{
    Matrix matrix = Matrix(3, 3);
    int* a = matrix[1];
    a[1] = 1;

    int* b = matrix[1];
    std::cout << b[1] << std::endl;

    return 0;
}