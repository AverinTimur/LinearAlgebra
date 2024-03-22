#include "HEADERS.h"
#include "program.cpp"
using namespace LinearAlgebra;

void print_matrix(Matrix matrix, int m, int n)
{
    std::cout << std::endl;
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

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

    try
    {
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
        print_matrix(matrix, 3, 3);

        std::cout << matrix.GetSubMatrixDeterminant(count, rows, coloms) << std::endl;
        std::cout << matrix.GetMinor(count, rows, coloms) << std::endl;
        std::cout << matrix.GetDeterminator() << std::endl;

        matrix.SwapRows(1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SwapRows(1, 2);

        matrix.SwapColomns(1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SwapColomns(1, 2);
        
        matrix.MultiplicatRow(4, 0);
        print_matrix(matrix, 3, 3);
        matrix.MultiplicatRow(0.25, 0);
        
        matrix.MultiplicatColomn(4, 0);
        print_matrix(matrix, 3, 3);
        matrix.MultiplicatColomn(0.25, 0);

        matrix.SumRow(0, 1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SumRow(0, 1, -2);
                
        matrix.SumColomn(0, 1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SumColomn(0, 1, -2);

    }
    catch(MatrixError error)
    {
        std::cout << std::endl;
        std::cout << "MATRIX ERROR" << std::endl;
        std::cout << error.name << std::endl;
        std::cout << error.text << std::endl;
    }

    return 0;
}