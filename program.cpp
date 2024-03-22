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
            float **array;

        public:
            // Base functions
            Matrix(int m, int n)
            {
                this->n = n;
                this->m = m;
                array = new float*;
                for(int i = 0; i < m; i++)
                {
                    array[i] = new float[n];
                    for(int j = 0; j < n; j++)
                    {
                        array[i][j] = 0;
                    }
                }
            }
            
            float* operator [](int i)
            {
                return array[i];
            }

            // Determinant calculation
            float GetSubMatrixDeterminant(int count, int rows[], int colomns[])
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
                float determinant;
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
            float GetMinor(int count, int rows[], int colomns[])
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
            float GetDeterminator()
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
            
            // Elementary operations
            void SwapRows(int first_row, int second_row)
            {
                float buffer[n];
                for(int i = 0; i < m; i++)
                {
                    buffer[i] = (*this)[first_row][i];
                }
                for(int i = 0; i < m; i++)
                {
                    (*this)[first_row][i] = (*this)[second_row][i];
                    (*this)[second_row][i] = buffer[i];
                }
            }
            void SwapColomns(int first_colomn, int second_colomn)
            {
                float buffer[n];
                for(int i = 0; i < m; i++)
                {
                    buffer[i] = (*this)[i][first_colomn];
                }
                for(int i = 0; i < m; i++)
                {
                    (*this)[i][first_colomn] = (*this)[i][second_colomn];
                    (*this)[i][second_colomn] = buffer[i];
                }
            }
            void MultiplicatRow(float k, int row)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[row][i] *= k;
                }
            }
            void MultiplicatColomn(float k, int colomn)
            {
                for(int i = 0; i < n; i++)
                {
                    (*this)[i][colomn] *= k;
                }
            }
            void SumRow(int old_row, int added_row, int k)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[old_row][i] += k * (*this)[added_row][i];
                }
            }
            void SumColomn(int old_colomn, int added_colomn, int k)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[i][old_colomn] += k * (*this)[i][added_colomn];
                }
            }
    };
}