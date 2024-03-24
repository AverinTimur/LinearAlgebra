#include "HEADERS.h"
#include "program.cpp"
using namespace LinearAlgebra;

void print_matrix(Matrix<float> matrix, int m, int n)
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
        Matrix matrix = Matrix<float>(3, 3);
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
        std::cout << matrix.GetDeterminant() << std::endl;

        matrix.SwapRows(1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SwapRows(1, 2);

        matrix.SwapColumns(1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SwapColumns(1, 2);
        
        matrix.MultiplyRow(4, 0);
        print_matrix(matrix, 3, 3);
        matrix.MultiplyRow(0.25, 0);
        
        matrix.MultiplyColumns(4, 0);
        print_matrix(matrix, 3, 3);
        matrix.MultiplyColumns(0.25, 0);

        matrix.SumRow(0, 1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SumRow(0, 1, -2);
                
        matrix.SumColumns(0, 1, 2);
        print_matrix(matrix, 3, 3);
        matrix.SumColumns(0, 1, -2);

        std::cout << matrix.GetRank() << std::endl;

        Matrix new_matrix = Matrix<float>(3, 2);
        k = 1;
        for(int i = 0; i < 2; i++)
        {
            for(int j = 0; j < 2; j++)
            {  
                new_matrix[i][j] = k;
                k++;
            }
        }
        new_matrix[0][1] = -2;
        new_matrix[2][0] = 1;
        new_matrix[2][1] = 1;

        float constant_terms[3]{0, 10, 3};
        if(new_matrix.IsSolutionExists(constant_terms))
        {
            std::cout << "True" << std::endl;
        }
        else
        {
            std::cout << "False" << std::endl;
        }
        
        constant_terms[2] = 2;
        if(new_matrix.IsSolutionExists(constant_terms))
        {
            std::cout << "True" << std::endl;
        }
        else
        {
            std::cout << "False" << std::endl;
        }

        new_matrix = Matrix<float>(3, 3);
        new_matrix[0][0] = 5;
        new_matrix[1][0] = 4;
        new_matrix[2][0] = 4;
        new_matrix[0][1] = 4;
        new_matrix[1][1] = 5;
        new_matrix[2][1] = 4;
        new_matrix[0][2] = 4;
        new_matrix[1][2] = 4;
        new_matrix[2][2] = 5;
        constant_terms[0] = 11;
        constant_terms[1] = 8;
        constant_terms[2] = 7;

        float* solution = new_matrix.CramerRule(constant_terms);
        for(int i = 0; i < 3; i++)
        {
            std::cout << solution[i] << " ";
        }
        std::cout << std::endl;
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