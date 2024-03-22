#include "HEADERS.h"

namespace LinearAlgebra
{
    class Matrix
    {
        private:
            int n, m; // count of coloms, rows
            int **array;

        public:
            Matrix(int n, int m)
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
    };
}