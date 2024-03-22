#include "HEADERS.h"
#include "program.cpp"
using namespace LinearAlgebra;

int main()
{
    int count;
    std::cin >> count;
    
    int rows[count], coloms[count];
    for(int i = 0; i < count; i++)
    {
        std::cin >> rows[i];
    }
    for(int i = 0; i < count; i++)
    {
        std::cin >> coloms[i];
    }

    Matrix matrix = Matrix(3, 3);
    int k = 0;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {  
            matrix[i][j] = k;
            k++;
        }
    }

    std::cout << matrix.GetSubMatrixDeterminant(count, rows, coloms) << std::endl;
    std::cout << matrix.GetMinor(count, rows, coloms) << std::endl;
    std::cout << matrix.GetDeterminator() << std::endl;

    return 0;
}