//
// Created by artoria on 5/13/20.
//

#ifndef MPI_COMMONTOOLS_H
#define MPI_COMMONTOOLS_H

#include <cstdint>
#include <fstream>

template<std::size_t n, std::size_t m>
constexpr void transpose( double array[n * m] ) {
    double new_array[m][n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            // Index in the original matrix.
            int index1 = i * m + j;

            // Index in the transpose matrix.
            // int index2 = j * m + i;
            new_array[j][i] = array[index1];
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            array[i * n + j] = new_array[i][j];
        }
    }
}

void write_out( double *T, int phi_size, int r_size, double dphi, double dr, const double R[2] );

#endif //MPI_COMMONTOOLS_H
