#include "HEADERS.h"

namespace LinearAlgebra
{
    class Matrix
    {
        private:
            int m, n; // count of rows, coloms
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
            int GetSubMatrixDeterminant(int count, int rows[], int coloms[])
            {
                // Error check
                for(int i = 0; i < count; i += 1)
                {
                    if(rows[i] >= m)
                    {
                        throw "ERROR";
                    }
                }
                for(int i = 0; i < count; i += 1)
                {
                    if(coloms[i] >= n)
                    {
                        throw "ERROR";
                    }
                }

                // Calculation
                int determinant;
                if(count == 1)
                {
                    determinant =  (*this)[rows[0]][coloms[0]];
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
                            determinant += (*this)[rows[i]][coloms[0]] * GetSubMatrixDeterminant((count - 1), new_rows, (coloms + 1));
                        }
                        else
                        {
                            determinant -= (*this)[rows[i]][coloms[0]] * GetSubMatrixDeterminant((count - 1), new_rows, (coloms + 1));
                        }
                    }
                }

                return determinant;
            }
            int GetMinor(int count, int rows[], int coloms[])
            {
                // Error check
                if(m != n)
                {
                    throw "ERROR";
                }

                int new_rows[n - count], new_coloms[n - count];
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
                    while(last_number < coloms[i])
                    {
                        new_coloms[j] = last_number;
                        last_number++;
                        j++;
                    }
                    last_number += 1;
                }
                while(last_number < n)
                {
                    new_coloms[j] = last_number;
                    last_number++;
                    j++;
                }

                return GetSubMatrixDeterminant(n - count, new_rows, new_coloms);
            }
            int GetDeterminator()
            {
                // Error check
                if(m != n)
                {
                    throw "ERROR";
                }

                int new_rows[n], new_coloms[n];
                for(int i = 0; i < n; i++)
                {
                    new_rows[i] = i;
                }
                for(int i = 0; i < n; i++)
                {
                    new_coloms[i] = i;
                }
                
                return GetSubMatrixDeterminant(n, new_rows, new_coloms);
            }
            
    };
}