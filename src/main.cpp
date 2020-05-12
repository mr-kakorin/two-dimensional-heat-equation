#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>


template<typename T, std::size_t N>
struct CompileTimeArray {
	constexpr explicit CompileTimeArray() : arr() {
	}

	[[nodiscard]] constexpr T get( int idx ) {
		return arr[idx];
	}

	[[nodiscard]] constexpr T get( int idx ) const {
		return arr[idx];
	}

	mutable T arr[N];
};

template<std::size_t num_r_steps, std::size_t num_phi_steps>
struct T0
		: public CompileTimeArray<double, num_r_steps * num_phi_steps> {
	constexpr T0() : CompileTimeArray<double, num_r_steps * num_phi_steps>() {
		init();
	}

	constexpr void init() {
		constexpr double phi1 = 0;
		constexpr double phi2 = M_PI * 2.;
		constexpr double Phi[2] = {phi1, phi2};
		constexpr double Tstart = 283.15;
		for (int i = 0; i < num_phi_steps; ++i) {
			for (int j = 0; j < num_r_steps; ++j) {
				if (j == 0) {
					CompileTimeArray<double, num_r_steps * num_phi_steps>::arr[num_phi_steps * i + j] = 293.15;
				} else if (j == num_r_steps - 1) {
					if (((j + 1) * Phi[1] / num_phi_steps >= 0 && (j + 1) * Phi[1] / num_phi_steps <= M_PI_2) ||
					    ((j + 1) * Phi[1] / num_phi_steps >= 3 * M_PI_2 &&
					     (j + 1) * Phi[1] / num_phi_steps <= 2 * M_PI)) {
						CompileTimeArray<double, num_r_steps * num_phi_steps>::arr[num_phi_steps * i + j] = 278.15;
					} else {
						CompileTimeArray<double, num_r_steps * num_phi_steps>::arr[num_phi_steps * i + j] = 283.15;
					}
				} else {
					CompileTimeArray<double, num_r_steps * num_phi_steps>::arr[num_phi_steps * i + j] = Tstart;
				}
			}
		}
	}

	constexpr double operator[]( int idx ) {
		return CompileTimeArray<double, num_r_steps * num_phi_steps>::get( idx );
	}

	constexpr double operator[]( int idx ) const {
		return CompileTimeArray<double, num_r_steps * num_phi_steps>::get( idx );
	}

	constexpr void swap( int idx, double x ) const {
		std::exchange( CompileTimeArray<double, num_r_steps * num_phi_steps>::arr[idx], x );
	}
};


void exchange( double *T0, const int &start, const int &end, const int &size, const int &layers,
               const int &rank, const int &num_pc ) {
	MPI_Allgather( T0 + start * size, size * layers, MPI_DOUBLE, T0, size * layers, MPI_DOUBLE, MPI_COMM_WORLD);
}

void end_exchange( const int &rank, const int &num_pc, double &maxerr ) {
	MPI_Status status;
	if (rank > 0) {
		MPI_Send((&maxerr), 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	} else {
		for (int i = 1; i < rank; ++i) {
			double tmp;
			MPI_Recv((&tmp), 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status );
			if (tmp > maxerr)
				maxerr = tmp;
		}
	}
}

template<std::size_t num_r_steps, std::size_t num_phi_steps>
void
first_stage( const double *T0, double *T1, const int &startPhi, const int &endPhi, const double &a, const double &dt,
             const double &dphi,
             const double &dr, const double R[2] ) {
	//coefficients
	double A[num_r_steps], B[num_r_steps];
	for (int j = startPhi; j < endPhi; ++j) {
		//double phi = (j + 1) * dphi;
		//start
		//T1[num_phi_steps * j + 0] = 283.15;//T0[num_phi_steps*j+0];
		//end
		//if ((phi >= 0 && phi <= M_PI_2) ||
		//    (phi >= 3 * M_PI_2 && phi <= 2 * M_PI)) {
		//	T1[num_phi_steps * j + num_r_steps - 1] = 278.15;
		//} else {
		//	T1[num_phi_steps * j + num_r_steps - 1] = 283.15;
		//}
		//T1[num_phi_steps*j + num_r_steps -1] = T0[num_phi_steps*j + num_r_steps -1];

		//coefficients beta and alpha
		A[0] = 0;
		B[0] = T0[0];
		double dr2 = dr * dr;
		double dta_dr2 = dt * a / dr2;
		for (int i = 0; i < num_r_steps - 1; ++i) {
			double ri = R[0] + (i + 1) * dr;
			double Ai = dta_dr2 * (1. - dr / (2. * ri));
			double Bi = -(1 + 2 * dta_dr2);
			double Ci = dta_dr2 * (1 + dr / (2. * ri));
			double Fi = -T0[i + 1];
			A[i + 1] = -Ci / (Bi + A[i] * Ai);
			B[i + 1] = (Fi - Ai * B[i]) / (Bi + A[i] * Ai);
		}
		for (int i = num_r_steps - 1; i >= 1; --i) {
			T1[num_phi_steps * j + i - 1] = A[i] * T1[num_phi_steps * j + i] + B[i];
		}
	}
}

void write_out( double *T, int phi_size, int r_size, double dphi, double dr, const double R[2] ) {
	std::ofstream out;
	out.open( "out.txt" );
	for (int i = 0; i < phi_size; ++i) {
		double phi = i * dphi;
		for (int j = 0; j < r_size; ++j) {
			double r = j * dr + R[0];
			out << r << " " << phi << " " << T[i * phi_size + j] << std::endl;
		}
		out << std::endl;
	}
	out.close();
}

template<std::size_t num_r_steps, std::size_t num_phi_steps>
void second_stage( const double *T0, double *T1, const int &startR, const int &endR, const double &a, const double &dt,
                   const double &dphi,
                   const double &dr, const double R[2] ) {
	double A[num_phi_steps], B[num_phi_steps], C[num_phi_steps];
	double P[num_phi_steps], Q[num_phi_steps];
	for (int j = startR; j < endR; ++j) {
		double ri = R[0] + (j + 1) * dr;
		double adt = a * dt;
		double ri2 = ri * ri;
		double dphi2 = dphi * dphi;
		double dphi2ri2_adt = ri2 * dphi2 / adt;

		double Ci = 2 + dphi2ri2_adt;
		A[0] = 1. / Ci;
		B[0] = dphi2ri2_adt * T0[j + 0 * num_r_steps];
		C[0] = 1. / Ci;
		for (int i = 1; i < num_phi_steps; ++i) {
			double Fi = dphi2ri2_adt * T0[j + i * num_r_steps];
			A[i] = 1. / (Ci - A[i - 1]);
			B[i] = (Fi + B[i - 1]) / (Ci - A[i - 1]);
			C[i] = C[i - 1] / (Ci - A[i - 1]);
		}

		P[num_phi_steps - 1] = B[num_phi_steps - 1];
		Q[num_phi_steps - 1] = A[num_phi_steps - 1] + C[num_phi_steps - 1];

		for (int i = num_phi_steps - 2; i >= 0; --i) {
			P[i] = A[i + 1] * P[i + 1] + B[i + 1];
			Q[i] = A[i + 1] * Q[i + 1] + C[i + 1];
		}

		T1[j + num_r_steps * (num_r_steps - 1)] = (B[num_phi_steps - 1] + A[num_phi_steps - 1] * P[0]) /
		                                          (1 - A[num_phi_steps - 1] * Q[0] - C[num_phi_steps - 1]);
		T1[j + 0 * num_r_steps] = T1[j + num_r_steps * (num_r_steps - 1)];
		for (int i = 1; i < num_phi_steps - 1; ++i) {
			T1[j + num_r_steps * i] = P[i] + T1[j + num_r_steps * (num_r_steps - 1)] * Q[i];
		}
	}
}

template<std::size_t n, std::size_t m>
constexpr void transpose( double array[n * m] ) {
	double new_array[m][n];
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			// Index in the original matrix.
			int index1 = i * n + j;

			// Index in the transpose matrix.
			// int index2 = j * m + i;
			new_array[j][i] = array[index1];
		}
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; ++j)
			array[i * m + j] = new_array[i][j];
	}
}

template<std::size_t num_r_steps, std::size_t num_phi_steps>
constexpr void init( double array[num_r_steps * num_phi_steps] ) {
	constexpr double phi1 = 0;
	constexpr double phi2 = M_PI * 2.;
	constexpr double Phi[2] = {phi1, phi2};
	constexpr double Tstart = 283.15;
	for (int i = 0; i < num_phi_steps; ++i) {
		for (int j = 0; j < num_r_steps; ++j) {
			if (j == 0) {
				array[num_phi_steps * i + j] = 293.15;
			} else if (j == num_r_steps - 1) {
				if (((j + 1) * Phi[1] / num_phi_steps >= 0 && (j + 1) * Phi[1] / num_phi_steps <= M_PI_2) ||
				    ((j + 1) * Phi[1] / num_phi_steps >= 3 * M_PI_2 &&
				     (j + 1) * Phi[1] / num_phi_steps <= 2 * M_PI)) {
					array[num_phi_steps * i + j] = 278.15;
				} else {
					array[num_phi_steps * i + j] = Tstart;
				}
			} else {
				array[num_phi_steps * i + j] = Tstart;
			}
		}
	}
}

template<int num_r_steps, int num_phi_steps>
constexpr void start_solve( const int &rank, const int &num_pc, const double R[2], const double Phi[2], const double a,
                            const double dt ) {
	double initial_temperature_matrix[num_r_steps * num_phi_steps];
	init<num_r_steps, num_phi_steps>( initial_temperature_matrix );
	float fTimeStart = 0;
	if (rank == 0)
		fTimeStart = clock() / (float) CLOCKS_PER_SEC;
	double dr = (R[1] - R[0]) / num_r_steps;

	int num_r_steps_rank = num_r_steps / num_pc;

	int startR = num_r_steps_rank * rank;
	int endR = startR + num_r_steps_rank;

	if (rank == num_pc - 1) {
		endR = num_r_steps;
		num_r_steps_rank = num_r_steps - startR;
	}

	const int r_size = num_r_steps_rank;

	int num_phi_steps_rank = num_phi_steps / num_pc;
	double dphi = (Phi[1]) / num_phi_steps;
	int startPhi = num_phi_steps_rank * rank;
	int endPhi = startPhi + num_phi_steps_rank;

	if (rank == num_pc - 1) {
		endPhi = num_phi_steps;
		num_phi_steps_rank = num_phi_steps - startPhi;
	}

	const int phi_size = num_phi_steps_rank;

	double T_previous_layer[num_phi_steps * num_r_steps];
	double T_current_layer[num_phi_steps * num_r_steps];

	for (int i = 0; i < num_phi_steps; ++i) {
		for (int j = 0; j < num_r_steps; ++j)
			T_previous_layer[num_phi_steps * i + j] = initial_temperature_matrix[num_phi_steps * i + j];
	}
	//delete T0;
	int l = 0;
	double nev = 0.;
	while (true) {
		first_stage<num_r_steps, num_phi_steps>( T_previous_layer, T_current_layer, startPhi, endPhi, a, dt, dphi,
		                                         dr,
		                                         R );
		std::swap( T_current_layer, T_previous_layer );
		MPI_Barrier(MPI_COMM_WORLD);
		exchange( T_previous_layer, startPhi, endPhi, num_r_steps, phi_size, rank, num_pc );

		second_stage<num_r_steps, num_phi_steps>( T_previous_layer, T_current_layer, startR, endR, a, dt, dphi,
		                                          dr, R );
		std::swap( T_current_layer, T_previous_layer );
		transpose<num_phi_steps, num_r_steps>( T_previous_layer );
		MPI_Barrier(MPI_COMM_WORLD);

		exchange( T_previous_layer, startR, endR, num_phi_steps, r_size, rank, num_pc );
		transpose<num_phi_steps, num_r_steps>( T_previous_layer );
		double err = 0.;
		for (int i = 0; i < num_phi_steps; ++i) {
			double terr = 0.;
			for (int j = 0; j < num_r_steps; ++j) {
				if (l > 0) {
					//double terr = std::abs(T0[num_phi_steps * i + j] - T_previous_layer[num_phi_steps * i + j]);
					if (std::abs( initial_temperature_matrix[num_phi_steps * i + j] -
					              T_previous_layer[num_phi_steps * i + j] ) > terr) {
						terr = std::abs( initial_temperature_matrix[num_phi_steps * i + j] -
						                 T_previous_layer[num_phi_steps * i + j] );
						//std::cout<<T0[num_phi_steps * i + j] <<" "<< T_previous_layer[num_phi_steps * i + j]<<std::endl;
					}
				}
				initial_temperature_matrix[num_phi_steps * i + j] = T_previous_layer[num_phi_steps * i + j];
			}
			if (terr > err)
				err = terr;
		}
		if (l % 1 == 0)
			std::cout << "Current nevyazka on process " << rank << " = " << err << std::endl;

		if (l++ > 1 && err < 1e-3) {
			nev = err;
			break;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	end_exchange( rank, num_pc, nev );

	if (rank == 0) {
		std::cout << "Time of static temperature distributio = " << (l + 1) * dt << std::endl;
		std::cout << "Static temperature distribution error = " << nev << std::endl;
		float fTimeStop = clock() / (float) CLOCKS_PER_SEC;
		std::cout << "Working time = " << -fTimeStart + fTimeStop << std::endl;
		write_out( T_previous_layer, num_phi_steps, num_r_steps, dphi, dr, R );
	}
}

//-np 4 /home/mrkakorin/CLionProjects/task2/cmake-build-debug/build/task2
int main( int argc, char **argv ) {
	//float fTimeStart = clock() / (float)CLOCKS_PER_SEC;
	constexpr double R1 = 0.5;
	constexpr double R2 = 2.9;
	constexpr double phi1 = 0;
	constexpr double phi2 = M_PI * 2.;
	constexpr double R[2] = {R1, R2};
	constexpr double Phi[2] = {phi1, phi2};


	constexpr int num_r_steps = 240;
	constexpr int num_phi_steps = 240;
	constexpr double dt{0.5};
	constexpr double a{5.33e-6};

	MPI_Init( &argc, &argv );
	int num_pc = 0, rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );
	MPI_Comm_size(MPI_COMM_WORLD, &num_pc );
	std::cout << num_pc << " " << rank << std::endl;
	//constexpr double T0[num_phi_steps * num_r_steps] = {};

	start_solve<num_r_steps, num_phi_steps>( rank, num_pc, R, Phi, a, dt );
	MPI_Finalize();
//  float fTimeStop = clock() / (float)CLOCKS_PER_SEC;
//  if(rank == 0)
//    std::cout << "Working time = " << -fTimeStart + fTimeStop << std::endl;
	//sin(pi*t)
	return 0;
}
