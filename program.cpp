#include "HEADERS.h"

namespace LinearAlgebra
{
    struct MatrixError
    {
        std::string name;
        std::string text;

        MatrixError(std::string name, std::string text)
        {
            this->name = name;
            this->text = text;
        }
    };
    class Matrix
    {
        private:
            int m, n; // count of rows, colomns
            int **array;

        public:
            // Base functions
            Matrix(int m, int n)
            {
                this->n = n;
                this->m = m;
                array = new int*;
                for(int i = 0; i < m; i++)
                {
                    array[i] = new int[n];
                    for(int j = 0; j < n; j++)
                    {
                        array[i][j] = 0;
                    }
                }
            }
            
            int* operator [](int i)
            {
                return array[i];
            }

            // Parameters calculation
            int GetSubMatrixDeterminant(int count, int rows[], int colomns[])
            {
                // Error check
                for(int i = 0; i < count; i += 1)
                {
                    if(rows[i] >= m)
                    {
                        throw MatrixError("Number out of Matrix", "Matrix don't contains row" + std::to_string(m));
                    }
                }
                for(int i = 0; i < count; i += 1)
                {
                    if(colomns[i] >= n)
                    {
                        throw MatrixError("Number out of Matrix", "Matrix don't contains colomn" + std::to_string(n));
                    }
                }

                // Calculation
                int determinant;
                if(count == 1)
                {
                    determinant =  (*this)[rows[0]][colomns[0]];
                }
                else
                {
                    determinant = 0;
                    for(int i = 0; i < count; i++)
                    {
                        int new_rows[count - 1];
                        for(int j = 0; j < count; j++)
                        {
                            if(j < i)
                            {
                                new_rows[j] = rows[j];
                            }
                            else
                            {
                                new_rows[j] = rows[j + 1];
                            }
                        }   

                        if(i % 2 == 0)
                        {
                            determinant += (*this)[rows[i]][colomns[0]] * GetSubMatrixDeterminant((count - 1), new_rows, (colomns + 1));
                        }
                        else
                        {
                            determinant -= (*this)[rows[i]][colomns[0]] * GetSubMatrixDeterminant((count - 1), new_rows, (colomns + 1));
                        }
                    }
                }

                return determinant;
            }
            int GetMinor(int count, int rows[], int colomns[])
            {
                // Error check
                if(m != n)
                {
                    throw MatrixError("Matrix isn't square", "Minor is determined only in square matrix");
                }

                int new_rows[n - count], new_colomns[n - count];
                int last_number = 0;
                int j = 0;
                for(int i = 0; i < count; i++)
                {
                    while(last_number < rows[i])
                    {
                        new_rows[j] = last_number;
                        last_number++;
                        j++;
                    }
                    last_number += 1;
                }
                while(last_number < n)
                {
                    new_rows[j] = last_number;
                    last_number++;
                    j++;
                }
                last_number = 0;
                j = 0;
                for(int i = 0; i < count; i++)
                {
                    while(last_number < colomns[i])
                    {
                        new_colomns[j] = last_number;
                        last_number++;
                        j++;
                    }
                    last_number += 1;
                }
                while(last_number < n)
                {
                    new_colomns[j] = last_number;
                    last_number++;
                    j++;
                }

                return GetSubMatrixDeterminant(n - count, new_rows, new_colomns);
            }
            int GetDeterminator()
            {
                // Error check
                if(m != n)
                {
                    throw MatrixError("Matrix isn't square", "Determinator is determined only in square matrix");
                }

                int new_rows[n], new_colomns[n];
                for(int i = 0; i < n; i++)
                {
                    new_rows[i] = i;
                }
                for(int i = 0; i < n; i++)
                {
                    new_colomns[i] = i;
                }
                
                return GetSubMatrixDeterminant(n, new_rows, new_colomns);
            }
            
    };
}