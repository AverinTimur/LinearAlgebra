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
            int m, n; // count of rows, columns
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
            Matrix Clone()
            {
                Matrix matrix = Matrix(m, n);

                for(int i = 0; i < m; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        matrix[i][j] = (*this)[i][j];
                    }
                }

                return matrix;
            }

            float* operator [](int i)
            {
                return array[i];
            }

            int GetRowOrder()
            {
                return m;
            }
            int GetColumnOrder()
            {
                return n;
            }

            // Determinant calculation
            float GetSubMatrixDeterminant(int count, int rows[], int columns[])
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
                    if(columns[i] >= n)
                    {
                        throw MatrixError("Number out of Matrix", "Matrix don't contains column" + std::to_string(n));
                    }
                }

                // Calculation
                std::sort(rows, rows + count, [](float a, float b){return a < b;});
                std::sort(columns, columns + count, [](float a, float b){return a < b;});

                float determinant;
                if(count == 1)
                {
                    determinant = (*this)[rows[0]][columns[0]];
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
                            determinant += (*this)[rows[i]][columns[0]] * GetSubMatrixDeterminant((count - 1), new_rows, (columns + 1));
                        }
                        else
                        {
                            determinant -= (*this)[rows[i]][columns[0]] * GetSubMatrixDeterminant((count - 1), new_rows, (columns + 1));
                        }
                    }
                }

                return determinant;
            }
            float GetMinor(int count, int rows[], int columns[])
            {
                // Error check
                if(m != n)
                {
                    throw MatrixError("Matrix isn't square", "Minor is determined only in square matrix");
                }

                int new_rows[n - count], new_columns[n - count];
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
                    while(last_number < columns[i])
                    {
                        new_columns[j] = last_number;
                        last_number++;
                        j++;
                    }
                    last_number += 1;
                }
                while(last_number < n)
                {
                    new_columns[j] = last_number;
                    last_number++;
                    j++;
                }

                return GetSubMatrixDeterminant(n - count, new_rows, new_columns);
            }
            float GetDeterminant()
            {
                // Error check
                if(m != n)
                {
                    throw MatrixError("Matrix isn't square", "Determinant is determined only in square matrix");
                }

                int new_rows[n], new_columns[n];
                for(int i = 0; i < n; i++)
                {
                    new_rows[i] = i;
                }
                for(int i = 0; i < n; i++)
                {
                    new_columns[i] = i;
                }
                
                return GetSubMatrixDeterminant(n, new_rows, new_columns);
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
            void SwapColumns(int first_column, int second_column)
            {
                float buffer[n];
                for(int i = 0; i < m; i++)
                {
                    buffer[i] = (*this)[i][first_column];
                }
                for(int i = 0; i < m; i++)
                {
                    (*this)[i][first_column] = (*this)[i][second_column];
                    (*this)[i][second_column] = buffer[i];
                }
            }
            void MultiplyRow(float k, int row)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[row][i] *= k;
                }
            }
            void MultiplyColumns(float k, int column)
            {
                for(int i = 0; i < n; i++)
                {
                    (*this)[i][column] *= k;
                }
            }
            void SumRow(int old_row, int added_row, int k)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[old_row][i] += k * (*this)[added_row][i];
                }
            }
            void SumColumns(int old_column, int added_column, int k)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[i][old_column] += k * (*this)[i][added_column];
                }
            }

            // Rank calculation
            int GetRank()
            {
                Matrix matrix = this->Clone();
                bool is_found;
                int rank = 0;

                do
                {
                    is_found = false;
                    for(int i = rank; i < m; i++)
                    {
                        if(is_found) break;
                        for(int j = rank; j < n; j++)
                        {
                            if(matrix[i][j] != 0)
                            {
                                is_found = true;
                                matrix.SwapRows(rank, i);
                                matrix.SwapColumns(rank, j);
                                break;
                            }
                        }
                    }

                    if(is_found)
                    {
                        matrix.MultiplyRow(1/matrix[rank][rank], rank);
                        for(int i = rank + 1; i < n; i++)
                        {
                            matrix.SumColumns(i, rank, -matrix[rank][i]);
                        }                        
                        for(int i = rank + 1; i < m; i++)
                        {
                            matrix.SumRow(i, rank, -matrix[i][rank]);
                        }

                        rank += 1;
                    }
                }
                while(is_found && rank < std::min(m, n));

                return rank;
            }

            // Linear system
            bool IsSolutionExists(float *constant_terms)
            {
                Matrix matrix = Matrix(m, n + 1);
                for(int i = 0; i < m; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        matrix[i][j] = (*this)[i][j];
                    }
                }
                for(int i = 0; i < m; i++)
                {
                    matrix[i][n] = constant_terms[i];
                }
                
                return this->GetRank() == matrix.GetRank();
            }

            int* GetBasis() // return r rows and r columns
            {
                int rank = GetRank();
                
                int *columns = new int[rank];
                int *rows = new int[2*rank];
                int count = 0;
                for(int j = 0; j < n; j++)
                {
                    columns[count] = j;
                    for(int i = rows[count - 1] + 1; i < m; i++)
                    {
                        rows[count] = i;
                        if(GetSubMatrixDeterminant(j + 1, rows, columns))
                        {
                            break;
                        }
                    }
                    count += 1;
                }

                for(int i = 0; i < rank; i++)
                {
                    rows[rank + i] = columns[i];
                }

                delete[] columns;
                return rows;
            }

            float* CramerRule(float *constant_terms)
            {
                int rank = GetRank();
                if(rank != n || !IsSolutionExists(constant_terms))
                {
                    throw MatrixError("Solution can't be find", "Solution isn't exists or it isn't unique");
                }

                int *basis = GetBasis();

                int determinant = GetSubMatrixDeterminant(rank, basis, (basis + rank));

                float *answer = new float[rank];
                float x;
                for(int j = 0; j < rank; j++)
                {
                    x = 0;
                    for(int i = 0; i < rank; i++)
                    { 
                        int new_rows[rank - 1];
                        int new_columns[rank - 1];
                        for(int t = 0; t < rank - 1; t++)
                        {
                            if(t < i)
                            {
                                new_rows[t] = basis[t];
                            }
                            else
                            {
                                new_rows[t] = basis[t + 1];
                            }  
                        }

                        for(int t = 0; t < rank - 1; t++)
                        {
                            if(t < j)
                            {
                                new_columns[t] = basis[rank + t];
                            }
                            else
                            {
                                new_columns[t] = basis[rank + t + 1];
                            }  
                        }
                        
                        if((basis[i] + basis[rank + j]) % 2 == 0)
                        {
                            x += constant_terms[basis[i]] * GetSubMatrixDeterminant(rank - 1, new_rows, new_columns);
                        }
                        else
                        {
                            x -= constant_terms[basis[i]] * GetSubMatrixDeterminant(rank - 1, new_rows, new_columns);
                        }
                    }
                    answer[j] = x / determinant;
                }

                delete[] basis;
                return answer;
            }   
    };
}