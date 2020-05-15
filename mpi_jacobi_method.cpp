#include <stdio.h>
#include <sys/time.h>
#include <cstdint>

#include <mpi.h>

#include <iostream>


const int32_t DIM = 10000;
const int32_t MAX_ITER = 10;
const double error_level = 1.0e-5;

double local_SE(double* first, double* second, int32_t start_idx, int32_t finish_idx) 
{
    double result = 0.0;
    for (int32_t i = start_idx; i < finish_idx; ++i) {
        result += (first[i] - second[i])*(first[i] - second[i]);
    }
    return result;
}

void mpi_init(int32_t cur_thread_index,
              int32_t threads_num,
              double** matrix,
              double* rhs,
              double* x_true,
              double* x,
              double* x_prev)
{
    int32_t start_idx = cur_thread_index * DIM / threads_num;
    int32_t finish_idx = (cur_thread_index + 1) * DIM / threads_num;
    int32_t part_size = finish_idx - start_idx;

    for(int32_t i = start_idx; i < finish_idx; ++i) {
        x_true[i] = i + 1.0;
        x[i] = 0.0;
        x_prev[i] = 0.0;
    }
    for (int32_t i = start_idx; i < finish_idx; ++i) {
        for (int32_t j = 0; j < DIM; ++j) {
            if (i == j) {
                matrix[i][j] = 100.0;
            } else if (i == j + 1 || i == j - 1) {
                matrix[i][j] = 1.0;
            } else {
                matrix[i][j] = 0.0;
            }
        }
    }
    // Sending part of x_true to other threads
    for (int32_t thread_idx = 0; thread_idx < threads_num; ++thread_idx) {
        if (cur_thread_index != thread_idx) {
            MPI_Send(&x_true[start_idx], part_size, MPI_DOUBLE, thread_idx, 0, MPI_COMM_WORLD);
        }
    }
    // Recieving other x_true parts from other threads
    for (int32_t thread_idx = 0; thread_idx < threads_num; ++thread_idx) {
        if (cur_thread_index != thread_idx) {
            int32_t tmp_start_idx = thread_idx * DIM / threads_num;
            int32_t tmp_finish_idx = (thread_idx + 1) * DIM / threads_num;
            int32_t tmp_part_size = tmp_finish_idx - tmp_start_idx;

            MPI_Recv(&x_true[tmp_start_idx], tmp_part_size, MPI_DOUBLE, thread_idx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    for (int32_t i = start_idx; i < finish_idx; ++i) {
        double tmp = 0.0;
        for (int32_t j = 0; j < DIM; ++j) {
            tmp += matrix[i][j] * x_true[j];
        }
        rhs[i] = tmp;
    }
}

void mpi_jacobi_method(int32_t cur_thread_index,
                       int32_t threads_num,
                       double** matrix,
                       double* rhs,
                       double* x,
                       double* x_prev)
{
    int32_t iter = 1;
    double error = 10.0;

    int32_t start_idx = cur_thread_index * DIM / threads_num;
    int32_t finish_idx = (cur_thread_index + 1) * DIM / threads_num;
    int32_t part_size = finish_idx - start_idx;    

    while ((error > error_level) && (iter < MAX_ITER)) {

        // Sending part of x_prev to other threads
        for (int32_t thread_idx = 0; thread_idx < threads_num; ++thread_idx) {
            if (cur_thread_index != thread_idx) {
                MPI_Send(&x_prev[start_idx], part_size, MPI_DOUBLE, thread_idx, 1, MPI_COMM_WORLD);
            }
        }
        // Recieving other x_prev parts from other threads
        for (int32_t thread_idx = 0; thread_idx < threads_num; ++thread_idx) {
            if (cur_thread_index != thread_idx) {
                int32_t tmp_start_idx = thread_idx * DIM / threads_num;
                int32_t tmp_finish_idx = (thread_idx + 1) * DIM / threads_num;
                int32_t tmp_part_size = tmp_finish_idx - tmp_start_idx;

                MPI_Recv(&x_prev[tmp_start_idx], tmp_part_size, MPI_DOUBLE, thread_idx, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        
        for (int32_t i = start_idx; i < finish_idx; ++i) {
            double sum_non_diag = 0.0;
            for (int32_t j = 0; j < DIM; ++j) {
                if (j != i) {
                    sum_non_diag += matrix[i][j] * x_prev[j];
                }
            }
            x[i] = (rhs[i] - sum_non_diag) / matrix[i][i];
        }
        double local_error = local_SE(x, x_prev, start_idx, finish_idx);
        MPI_Allreduce(&local_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        error = error / DIM;

        iter++;
        for (int32_t i = start_idx; i < finish_idx; ++i) {
            x_prev[i] = x[i];
        }
    }
}

int main()
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int32_t threads_num;
    MPI_Comm_size(MPI_COMM_WORLD, &threads_num);

    if (threads_num <= 0) {
        MPI_Finalize();
        std::cout << "Error is occured: number of processors is non-positive" << std::endl;
        return 0;
    }

    // Get the index of current process
    int32_t cur_thread_index;
    MPI_Comm_rank(MPI_COMM_WORLD, &cur_thread_index);    

    int32_t start_idx = cur_thread_index * DIM / threads_num;
    int32_t finish_idx = (cur_thread_index + 1) * DIM / threads_num;
    int32_t part_size = finish_idx - start_idx;

    double** matrix = new double*[DIM];
    for (int32_t i = 0; i < DIM; ++i) {
        matrix[i] = new double[DIM];
    }
    double* rhs = new double[DIM];

    double* x_true = new double[DIM];
    double* x = new double[DIM];
    double* x_prev = new double[DIM];

    mpi_init(cur_thread_index, threads_num, matrix, rhs, x_true, x, x_prev);

    if (cur_thread_index == 0) {    
        printf("Calculation is running at the %d thread(s)\n", threads_num);
    }

    struct timeval tv;

    gettimeofday(&tv, NULL);
    time_t start_sec = tv.tv_sec;
    suseconds_t start_usec = tv.tv_usec;

    mpi_jacobi_method(cur_thread_index, threads_num, matrix, rhs, x, x_prev);
    
    gettimeofday(&tv, NULL);
    time_t finish_sec = tv.tv_sec;
    suseconds_t finish_usec = tv.tv_usec;    

    // Sending part of x to other threads
    for (int32_t thread_idx = 0; thread_idx < threads_num; ++thread_idx) {
        if (cur_thread_index != thread_idx) {
            MPI_Send(&x[start_idx], part_size, MPI_DOUBLE, thread_idx, 2, MPI_COMM_WORLD);
        }
    }
    // Recieving other x parts from other threads
    for (int32_t thread_idx = 0; thread_idx < threads_num; ++thread_idx) {
        if (cur_thread_index != thread_idx) {
            int32_t tmp_start_idx = thread_idx * DIM / threads_num;
            int32_t tmp_finish_idx = (thread_idx + 1) * DIM / threads_num;
            int32_t tmp_part_size = tmp_finish_idx - tmp_start_idx;

            MPI_Recv(&x[tmp_start_idx], tmp_part_size, MPI_DOUBLE, thread_idx, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    /*if (cur_thread_index == 0) {
        printf("x vector is [");
        for (int32_t i = 0; i < DIM; ++i) {
            printf(" %f", x[i]);
        }
        printf(" ]\n");

        printf("x_true vector is [");
        for (int32_t i = 0; i < DIM; ++i) {
            printf(" %f", x_true[i]);
        }
        printf(" ]\n");
    }*/

    double final_error = 0.0;
    double local_final_error = local_SE(x, x_true, start_idx, finish_idx);
    MPI_Allreduce(&local_final_error, &final_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    final_error = final_error / DIM;
    if (cur_thread_index == 0) {
        printf("Final mean squared error is %e\n", final_error);
        printf("Time of calculation is %f sec\n", (finish_sec - start_sec) + (finish_usec - start_usec)*1.0e-6);
    }

    for (int32_t i = 0; i < DIM; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] rhs;

    delete[] x_true;
    delete[] x;
    delete[] x_prev;

    MPI_Finalize();

    return 0;
}
