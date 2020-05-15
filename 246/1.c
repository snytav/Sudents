#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

#define process 2
#define N 5520

struct timeval tv1, tv2;

int main() {

	gettimeofday(&tv1, NULL);
	double tau = 0.00001;
	const double epsilon = 0.0000001;
	double A[N * N];
	double x_n[N];
	double b[N];
	double x_n_1[N];

	for (long i = 0; i < N; i++) {
		for (long j = 0; j < N; j++) {
			A[i * N + j] = 1.0;
			if (i == j) {
				A[i * N + j] = 2.0;
			}
		}
	}

	for (long i = 0; i < N; i++) {
		b[i] = N + 1;
	}

	for (long i = 0; i < N; i++) {
		x_n[i] = 0;
	}


	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	double norm_Axn_b = 0;
	double norm_b = 0;
	
	MPI_Status status;

	for (long i = 0; i < N; i++) {
		norm_b += b[i] * b[i];
	}


	norm_b = sqrt(norm_b);

	char converged = 0;
	long itter = 0;

	int first = rank * (int)(N / process);
	int back = (rank + 1) * (int)(N / process);

	while (!converged) {

		double sub_norm_Axn_b = 0;

		if (rank == 0) {



				for (long i = first; i < back; i++) {
					double sub = 0;
					for (long j = 0; j < N; j++) {
						sub += A[i * N + j] * x_n[j];
					}
					x_n_1[i] = x_n[i] - tau * (sub - b[i]);
				}

				for (int k = 0; k < process; k++) {
					if (k != rank) {
						MPI_Send(&x_n_1[first], back - first, MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
						MPI_Recv(&x_n_1[k * (int)N / process], ((k+1) * (int)N) / process - (k * (int)N / process), MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
					}
				}

				double sub_norm_Axn_b = 0;
				memcpy(x_n, x_n_1, N * sizeof(double));

				for (long i = 0; i < N; i++) {
					double sub = 0;
					for (long j = 0; j < N; j++) {
						sub += (A[i * N + j]) * x_n[j];
					}
					sub_norm_Axn_b += (sub - b[i]) * (sub - b[i]);
				}


			double norm_from_another;

			for (int k = 0; k < process; k++) {
				if (k != rank) {
					MPI_Recv(&norm_from_another, 1, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
					sub_norm_Axn_b += norm_from_another;
				}
			}

			norm_Axn_b = sqrt(sub_norm_Axn_b);

			if (norm_Axn_b / norm_b < epsilon) {
				converged = 1;
			}

			for (int k = 0; k < process; k++) {
				MPI_Send(&converged, 1, MPI_INT, k, 1, MPI_COMM_WORLD);
			}

		}

		for (int r = 1; r < process; r++) {
			if (r == rank) {
				for (int i = first; i < back; i++) {
					
					double sub = 0;
					for (int j = 0; j < N; j++) {
						sub += *(A + i * N + j) * x_n[j];
					}
					x_n1[i] = x_n[i] - tau * (sub - b[i]);
				}

				for (int k = 0; k < process; k++) {
					if (k != rank) {
						MPI_Send(&x_n_1[first], back - first, MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
						MPI_Recv(&x_n_1[k * (int)N / process], ((k + 1) * (int)N) / process - (k * (int)N / process), MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
					}
				}

				memcpy(x_n, x_n_1, N * sizeof(double));

				for (long i = 0; i < N; i++) {
					double sub = 0;
					for (long j = 0; j < N; j++) {
						sub += (A[i * N + j]) * x_n[j];
					}
					sub_norm_Axn_b += (sub - b[i]) * (sub - b[i]);
				}

				MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
				MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

			}

		}

		itter++;

	}


	MPI_Finalize();

	printf("Itters: %i\n", itter);

	gettimeofday(&tv2, NULL);

	double dt_sec = (tv2.tv_sec - tv1.tv_sec);
	double dt_usec = (tv2.tv_usec - tv1.tv_usec);
	double dt = dt_sec + 1e-6 * dt_usec;
	printf("time diff %e \n", dt);




	return 0;
}

