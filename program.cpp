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

    template<typename T>
    class Matrix
    {
        private:
            int m, n; // count of rows, columns
            T **array;

        public:
            // Base functions
            Matrix(int m, int n)
            {
                this->n = n;
                this->m = m;
                array = new T*;
                for(int i = 0; i < m; i++)
                {
                    array[i] = new T[n];
                    for(int j = 0; j < n; j++)
                    {
                        array[i][j] = 0;
                    }
                }
            }
            void Dispose()
            {
                for(int i = 0; i < n; i++)
                {
                    delete[] array[i];
                }
                delete[] array;
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

            T* operator [](int i)
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
            T GetSubMatrixDeterminant(int count, int rows[], int columns[])
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
                std::sort(rows, rows + count, [](T a, T b){return a < b;});
                std::sort(columns, columns + count, [](T a, T b){return a < b;});

                T determinant;
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
            T GetMinor(int count, int rows[], int columns[])
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
            T GetDeterminant()
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
                if(first_row != second_row)
                {
                    T buffer[n];
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
            }
            void SwapColumns(int first_column, int second_column)
            {
                if(first_column != second_column)
                {
                    T buffer[n];
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
            }
            void MultiplyRow(T k, int row)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[row][i] *= k;
                }
            }
            void MultiplyColumns(T k, int column)
            {
                for(int i = 0; i < n; i++)
                {
                    (*this)[i][column] *= k;
                }
            }
            void SumRow(int old_row, int added_row, T k)
            {
                for(int i = 0; i < m; i++)
                {
                    (*this)[old_row][i] += k * (*this)[added_row][i];
                }
            }
            void SumColumns(int old_column, int added_column, T k)
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
            bool IsSolutionExists(T *constant_terms)
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
                
                int columns[rank];
                int *rows = new int[2*rank]{0};
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

                //delete[] columns;
                return rows;
            }

            T* CramerRule(T *constant_terms)
            {
                int rank = GetRank();
                if(rank != n || !IsSolutionExists(constant_terms))
                {
                    throw MatrixError("Solution can't be find", "Solution isn't exists or it isn't unique");
                }

                int *basis = GetBasis();

                int determinant = GetSubMatrixDeterminant(rank, basis, (basis + rank));

                T *answer = new T[rank];
                T x;
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

            T* GetTriangleAnalog(T *constant_terms)
            {
                bool is_string_null;
                T buffer;
                T *constant_terms_clone = new T[m];
                for(int i = 0; i < m; i++)
                {
                    constant_terms_clone[i] = constant_terms[i];
                }

                for(int i = 0; i < n - 1; i++)
                {
                    is_string_null = true;
                    for(int k = i; k < n; k++)
                    {
                        if((*this)[k][i] != 0)
                        {
                            is_string_null = false;
                            SwapRows(k, i);
                            buffer = constant_terms_clone[k];
                            constant_terms_clone[k] = constant_terms_clone[i];
                            constant_terms_clone[i] = buffer;
                            break;
                        }
                    }
                    
                    if(!is_string_null)
                    {
                        for(int j = i + 1; j < m; j++)
                        {
                            constant_terms_clone[j] -= (*this)[j][i]/(*this)[i][i]*constant_terms_clone[i];
                            SumRow(j, i, -(*this)[j][i]/(*this)[i][i]);
                            int a = 0;
                        }
                    }
                }

                return constant_terms_clone;
            }

            T* GaussianElimination(T *constant_terms)
            {
                T*constant_terms_clone = GetTriangleAnalog(constant_terms);
                
                bool is_string_null = true;
                for(int i = 0; i < m; i++)
                {
                    if((*this)[n - 1][i] != 0)
                    {
                        is_string_null = false;
                        break;
                    }
                }

                if(is_string_null)
                {
                    throw MatrixError("It isn't diagonal form", "Diagonal form because it's null string on matrix");
                }

                for(int i = n - 1; i >= 0; i--)
                {
                    int a = 0;
                    for(int j = 0; j < i; j++)
                    {
                        if((*this)[j][i] != 0)
                        {
                            constant_terms_clone[j] -= ((*this)[j][i]/(*this)[i][i])*constant_terms_clone[i];
                            SumRow(j, i, -(*this)[j][i]/(*this)[i][i]);
                        }
                    }
                }

                return constant_terms_clone;
            }

            T* GetSolutionByGaussianElimination(T *constant_terms)
            {
                Matrix clone_matrix = Clone();
                //T *constant_terms_clone = clone_matrix.GetTriangleAnalog(constant_terms);
                T *constant_terms_clone = clone_matrix.GaussianElimination(constant_terms);
                T *solution = new T[n];

                for(int i = 0; i < n; i++)
                {
                    if((*this)[i][i] == 0)
                    {
                        throw MatrixError("Solution can't be find", "Solution isn't exists or it isn't unique");
                    }
                    else
                    {
                        solution[i] = constant_terms_clone[i]/(clone_matrix[i][i]);
                    }
                }

                return solution;
            }
    };
}