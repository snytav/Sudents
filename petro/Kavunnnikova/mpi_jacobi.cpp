#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>
 
 
int main () {

    const int task_size = 2000;
    const double epsilon = 1.0e-5;

    double sum_non_diag = 0.0;
    double error = 50.0;
    int iter = 1;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Status status;  

    int cur_start = world_rank * task_size / world_size;
    int cur_finish = (world_rank + 1) * task_size / world_size;
    int cur_size = task_size / world_size;
    double cur_error = 0.0;


    /*       Initialization         */
    double** matrix = new double*[task_size];
    for (int i = 0; i < task_size; ++i) {
        matrix[i] = new double[task_size];
    }
    double* rhs = new double[task_size];
    double* x_exact = new double[task_size];
    double* x = new double[task_size];
    double* x_old = new double[task_size];

    for(int i = cur_start; i < cur_finish; ++i) {
        x_exact[i] = i;
        x[i] = 0.0;
        x_old[i] = 0.0;
    }
    for (int i = cur_start; i < cur_finish; ++i) {
        for (int j = 0; j < task_size; ++j) {
            if (i == j) {
                matrix[i][j] = 30.0;
            } else if (i == j - 1 || i == j + 1) {
                matrix[i][j] = 10.0;
            } else {
                matrix[i][j] = 0.0;
            }
        }
    }
    for (int i = 0; i < world_size; ++i) {
        if (i != world_rank) {
            MPI_Send(&x_exact[cur_start], cur_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }
    for (int i = 0; i < world_size; ++i) {
        if (i != world_rank) {
            MPI_Recv(&x_exact[i*cur_size], cur_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
        }
    }
    for (int i = cur_start; i < cur_finish; ++i) {
        for (int j = 0; j < task_size; ++j) {
            rhs[i] += matrix[i][j] * x_exact[j];
        }
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    time_t start_sec = tv.tv_sec;
    suseconds_t start_usec = tv.tv_usec;

    /*       Calculation         */

    while (error > epsilon) {
        for (int i = 0; i < world_size; ++i) {
            if (i != world_rank) {
                MPI_Send(&x_old[cur_start], cur_size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            }
        }
        for (int i = 0; i < world_size; ++i) {
            if (i != world_rank) {
                MPI_Recv(&x_old[i*cur_size], cur_size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
            }
        }
        for (int i = cur_start; i < cur_finish; ++i) {
            sum_non_diag = 0.0;
            for (int j = 0; j < task_size; ++j) {
                if (j != i) {
                    sum_non_diag = sum_non_diag + matrix[i][j]*x_old[j];
                }
            }
            x[i] = (rhs[i] - sum_non_diag) / matrix[i][i];
        }

        cur_error = 0.0;
        for (int i = cur_start; i < cur_finish; ++i) {
            cur_error = cur_error + (x[i] - x_old[i])*(x[i] - x_old[i]);
        }
        cur_error = sqrt(cur_error);
        MPI_Allreduce(&cur_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        error = error / task_size;

        for (int i = cur_start; i < cur_finish; ++i) {
            x_old[i] = x[i];
        }

        iter = iter + 1;
    }

    gettimeofday(&tv, NULL);
    time_t finish_sec = tv.tv_sec;
    suseconds_t finish_usec = tv.tv_usec;

    if (world_rank == 0) {
        printf("Time: %f\n", (finish_sec - start_sec) + (finish_usec - start_usec)*1.0e-6);
        printf("Error: %e\n", error);
        printf("Number of iterations: %d\n", iter);
    }
    
    for (int i = 0; i < task_size; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] rhs;

    delete[] x_exact;
    delete[] x;
    delete[] x_old;

    MPI_Finalize();

    return 0;
}