#include <fstream>
#include <math.h>
#include <cstdint>
#include <sys/time.h>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>

#include <mpi.h>

/* Calcualtion parameters */

const int32_t time_steps = 10;
const int32_t x_steps = 1000;
const int32_t y_steps = 500;

bool DEBUG = false;
bool WRITE_RESULT = true;

/* Physical parameters */

double min_x = 0.0;                                         // [m]
double max_x = 1.0;                                         // [m]
double min_y = 0.0;                                         // [m]
double max_y = 1.0;                                         // [m]

double d_x = (max_x - min_x) / x_steps;
double d_y = (max_y - min_y) / y_steps;

double alpha = 1.0;                                         // [m / s^2]
double start_time = 0.0;                                    // [s]

double final_time = std::min(d_x*d_x, d_y*d_y) / alpha;     // [s]

double d_t = (final_time - start_time) / time_steps;

double T_ampl = 1.0;                                        // [K]
double left_temperature = 0.0;                              // [K]
double right_temperature = 0.0;                             // [K]
double bottom_temperature = 0.0;                            // [K]
double top_temperature = 0.0;                               // [K]

/* Additional constants */

const double error_level = 1e-10;
const int32_t max_iter = 100;
const double PI = 3.1415926535;

const int32_t x_points = x_steps + 1;
const int32_t y_points = y_steps + 1;
const int32_t time_points = time_steps + 1;

int32_t index(int32_t i, int32_t j) {
    return i * y_points + j;
}

double exact_temperature(double x, double y, double t) {
    return T_ampl * sin(PI*x) * sin(PI*y) * exp(-1.0*t);
}

double F(double x, double y, double t) {
    return T_ampl * (2*alpha*PI*PI - 1) * sin(PI*x) * sin(PI*y) * exp(-1.0*t);
}

void init_coordinate_mesh(double* x, double* y) {
    for (int32_t i = 0; i < x_points; ++i) {
        x[i] = min_x + d_x * i;
    }
    for (int32_t i = 0; i < y_points; ++i) {
        y[i] = min_y + d_y * i;
    }
}

void update_exact_temperature(double* exact_temp,
                              double* x,
                              double* y,
                              double current_time,
                              int32_t cur_proc,
                              int32_t proc_num)
{
    int32_t start_i = cur_proc * x_points / proc_num;
    int32_t finish_i = (cur_proc + 1) * x_points / proc_num;
    int32_t part_size = finish_i - start_i;

    for (int32_t i = start_i; i < finish_i; ++i) {
        for (int32_t j = 0; j < y_points; ++j) {
            exact_temp[index(i, j)] = exact_temperature(x[i], y[j], current_time);
        }
    }
}

void init_temperature(double* temperature,
                      double* x,
                      double* y,
                      int32_t cur_proc,
                      int32_t proc_num)
{
    int32_t start_i = cur_proc * x_points / proc_num;
    int32_t finish_i = (cur_proc + 1) * x_points / proc_num;
    int32_t part_size = finish_i - start_i;

    for (int32_t i = start_i; i < finish_i; ++i) {
        for (int32_t j = 0; j < y_points; ++j) {
            temperature[index(i, j)] = exact_temperature(x[i], y[j], 0.0);
        }
    }
}

bool is_border(int32_t i, int32_t j) {
    return (i == 0 || j == 0 || i == x_points - 1 || j == y_points - 1);
}

void apply_boundary_condition(double* temperature,
                              double* x,
                              double* y,
                              double current_time,
                              int32_t cur_proc,
                              int32_t proc_num)
{
    int32_t start_i = cur_proc * x_points / proc_num;
    int32_t finish_i = (cur_proc + 1) * x_points / proc_num;
    int32_t part_size = finish_i - start_i;

    for (int32_t i = start_i; i < finish_i; ++i) {
        for (int32_t j = 0; j < y_points; ++j) {
            if (is_border(i, j)) {
                temperature[index(i, j)] = exact_temperature(x[i], y[j], current_time);
            }
        }
    }
}

void calculate_temperature(double* temperature,
                           double* x,
                           double* y,
                           double current_time,
                           int32_t cur_proc,
                           int32_t proc_num)
{   
    int32_t start_i = cur_proc * x_points / proc_num;
    int32_t finish_i = (cur_proc + 1) * x_points / proc_num;
    int32_t part_size = finish_i - start_i;

    double* new_temperature = new double [x_points * y_points];

    if (DEBUG && cur_proc == 0) {
        printf("New temperature array is initialized\n");
    }

    double* to_next_proc = new double [y_points];
    double* to_prev_proc = new double [y_points];

    for (int32_t j = 0; j < y_points; ++j) {
        to_next_proc[j] = temperature[index(finish_i - 1, j)];
        to_prev_proc[j] = temperature[index(start_i, j)];
    }
    if (DEBUG && cur_proc == 0) {
        printf("Temporary arrays are filled\n");
    }

    // Sending part of tempearture to next process (if exists)
    if (cur_proc != proc_num - 1) {
        MPI_Send(to_next_proc, y_points, MPI_DOUBLE, cur_proc + 1, 1, MPI_COMM_WORLD);
        if (DEBUG) {
            printf("Process %d sent data to process %d\n", cur_proc, cur_proc + 1);
        }
    }

    // Sending part of tempearture to prev process (if exists)
    if (cur_proc != 0) {
        MPI_Send(to_prev_proc, y_points, MPI_DOUBLE, cur_proc - 1, 0, MPI_COMM_WORLD);
        if (DEBUG) {
            printf("Process %d sent data to process %d\n", cur_proc, cur_proc - 1);
        }
    }

    double* from_next_proc = new double [y_points];
    double* from_prev_proc = new double [y_points];

    MPI_Status status;
    // Recieving part of temperature from next process (if exists)
    if (cur_proc != proc_num - 1) {
        MPI_Recv(from_next_proc, y_points, MPI_DOUBLE, cur_proc + 1, 0, MPI_COMM_WORLD, &status);
        if (DEBUG) {
            printf("Process %d recieved data from process %d\n", cur_proc, cur_proc + 1);
        }
    }

    // Recieving part of temperature from prev process (if exists)
    if (cur_proc != 0) {
        MPI_Recv(from_prev_proc, y_points, MPI_DOUBLE, cur_proc - 1, 1, MPI_COMM_WORLD, &status);
        if (DEBUG) {
            printf("Process %d recieved data from process %d\n", cur_proc, cur_proc - 1);
        }
    }

    apply_boundary_condition(new_temperature, x, y, current_time, cur_proc, proc_num);
    if (DEBUG && cur_proc == 0) {
        printf("Boundary conditions are applied\n");
    }
    
    for (int32_t i = start_i + 1; i < finish_i - 1; ++i) {
        for (int32_t j = 1; j < y_points - 1; ++j) {
            new_temperature[index(i, j)] = (1.0 - 2.0*d_t*alpha*(1.0 / (d_x*d_x) + 1.0 / (d_y*d_y))) * temperature[index(i, j)]
                                            + d_t*alpha / (d_x*d_x) * (temperature[index(i + 1, j)] + temperature[index(i - 1, j)])
                                            + d_t*alpha / (d_y*d_y) * (temperature[index(i, j + 1)] + temperature[index(i, j - 1)])
                                            + d_t*F(x[i], y[j], current_time);
        }
    }

    if (cur_proc != 0) {
        int32_t i = start_i;
        for (int32_t j = 1; j < y_points - 1; ++j) {
            new_temperature[index(i, j)] = (1.0 - 2.0*d_t*alpha*(1.0 / (d_x*d_x) + 1.0 / (d_y*d_y))) * temperature[index(i, j)]
                                            + d_t*alpha / (d_x*d_x) * (temperature[index(i + 1, j)] + from_prev_proc[j])
                                            + d_t*alpha / (d_y*d_y) * (temperature[index(i, j + 1)] + temperature[index(i, j - 1)])
                                            + d_t*F(x[i], y[j], current_time);
        }
    }

    if(cur_proc != proc_num - 1) {
        int32_t i = finish_i - 1;
        for (int32_t j = 1; j < y_points - 1; ++j) {
            new_temperature[index(i, j)] = (1.0 - 2.0*d_t*alpha*(1.0 / (d_x*d_x) + 1.0 / (d_y*d_y))) * temperature[index(i, j)]
                                            + d_t*alpha / (d_x*d_x) * (from_next_proc[j] + temperature[index(i - 1, j)])
                                            + d_t*alpha / (d_y*d_y) * (temperature[index(i, j + 1)] + temperature[index(i, j - 1)])
                                            + d_t*F(x[i], y[j], current_time);
        }
    }

    if (DEBUG && cur_proc == 0) {
        printf("New temperature is calculated\n");
    }

    for (int32_t i = start_i; i < finish_i; ++i) {
        for (int32_t j = 0; j < y_points; ++j) {
            temperature[index(i, j)] = new_temperature[index(i, j)];
        }
    }
    if (DEBUG && cur_proc == 0) {
        printf("Temperature value is updated\n");
    }

    delete[] to_next_proc;
    delete[] to_prev_proc;
    delete[] from_prev_proc;
    delete[] from_next_proc;
    
    delete[] new_temperature;

    if (DEBUG && cur_proc == 0) {
        printf("New temperature array is deleted\n");
    }
}

double SE(double* first, double* second,
          int32_t i_start, int32_t i_finish,
          int32_t j_start, int32_t j_finish) {
    double result = 0.0;
    for (int32_t i = i_start; i < i_finish; ++i) {
        for (int32_t j = j_start; j < j_finish; ++j) {
            result += (first[index(i, j)] - second[index(i, j)])*(first[index(i, j)] - second[index(i, j)]);
        }
    }
    return result;
}

void write_to_file(std::string filename,
                   double* temperature,
                   double* x, double* y) {
    FILE* file = fopen(filename.c_str(), "w");

    fprintf(file, "%16s", "x\\y");
    for (int32_t j = 0; j < y_points; ++j) {
        fprintf(file, "%16.3f", y[j]);
    }
    fprintf(file, "\n");
    for (int32_t i = 0; i < x_points; ++i) {
        fprintf(file, "%16.3f", x[i]);
        for (int32_t j = 0; j < y_points; ++j) {
            fprintf(file, "%16.3e", temperature[index(i, j)]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

int main(int argc, char** argv)
{
    struct timeval start_tv;
    gettimeofday(&start_tv, NULL);

    // Init phase

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int32_t proc_num;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

    if (proc_num <= 0) {
        MPI_Finalize();
	std::cout << "Error is occured: number of processors is non-positive" << std::endl;
        return 0;
    }

    // Get the index of current process
    int32_t cur_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &cur_proc);

    int32_t start_i = cur_proc * x_points / proc_num;
    int32_t finish_i = (cur_proc + 1) * x_points / proc_num;
    int32_t part_size = finish_i - start_i;    

    double* x = new double [x_points];
    double* y = new double [y_points];
    init_coordinate_mesh(x, y);

    if (DEBUG && cur_proc == 0) {
        printf("x array is initialized (%.3f kB)\n", x_points * sizeof(double) / 1024.0 );
        printf("y array is initialized (%.3f kB)\n", y_points * sizeof(double) / 1024.0 );
    }

    double* temperature = new double [x_points * y_points];
    double* exact_temperature = new double [x_points * y_points];

    init_temperature(temperature, x, y, cur_proc, proc_num);
    update_exact_temperature(exact_temperature, x, y, 0.0, cur_proc, proc_num);

    if (DEBUG && cur_proc == 0) {
        printf("Temperaure array is initialized (%.3f MB)\n", x_points * y_points * sizeof(double) / (1024.0 * 1024.0));
    }

    if (WRITE_RESULT && cur_proc == 0) {
        system("rm -r results");
        system("mkdir results");
        system("mkdir results/exact_temperature");
        system("mkdir results/temperature");
        write_to_file(std::string("results/temperature/initial.txt"), temperature, x, y);
        write_to_file(std::string("results/exact_temperature/initial.txt"), exact_temperature, x, y);
    } 

    struct timeval init_tv;
    gettimeofday(&init_tv, NULL);

    double init_time = (init_tv.tv_sec - start_tv.tv_sec) + (init_tv.tv_usec - start_tv.tv_usec)*1.0e-6;

    double solution_time = 0.0;
    double write_time = 0.0;

    // Solving phase

    for (int32_t time_step = 1; time_step <= time_steps; ++time_step) {

        struct timeval ref_tv;
        gettimeofday(&ref_tv, NULL);

        double current_time = time_step * d_t;
        
        if (cur_proc == 0) {
            printf("\nCurrent time = %.3e sec\n", current_time);
        }

        calculate_temperature(temperature, x, y, current_time, cur_proc, proc_num);
        if (DEBUG && cur_proc == 0) {
            printf("Temperature is calculated\n");
        }
        update_exact_temperature(exact_temperature, x, y, current_time, cur_proc, proc_num);
        if (DEBUG && cur_proc == 0) {
            printf("Exact temperature is updated\n");
        }

        double final_error = 0.0;
        double local_final_error = SE(temperature, exact_temperature, start_i, finish_i, 0, y_points);
        if (DEBUG) {
            printf("Local error with exact solution at the %d process = %.3e\n", cur_proc, local_final_error);
        }
        MPI_Allreduce(&local_final_error, &final_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
        final_error = final_error / (x_points * y_points);
        if (cur_proc == 0) {
            printf("Error with exact solution = %.3e\n", final_error);
        }

        struct timeval sol_tv;
        gettimeofday(&sol_tv, NULL);

        if (WRITE_RESULT) {

            if (cur_proc != 0) {
                MPI_Send(&temperature[index(start_i, 0)], part_size * y_points, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                MPI_Send(&exact_temperature[index(start_i, 0)], part_size * y_points, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
            }

            if (cur_proc == 0) {
                for (int32_t proc_idx = 1; proc_idx < proc_num; ++proc_idx) {

                    int32_t tmp_start_i = proc_idx * x_points / proc_num;
                    int32_t tmp_finish_i = (proc_idx + 1) * x_points / proc_num;
                    int32_t tmp_part_size = tmp_finish_i - tmp_start_i;

                    MPI_Recv(&temperature[index(tmp_start_i, 0)], tmp_part_size * y_points, MPI_DOUBLE, proc_idx, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&exact_temperature[index(tmp_start_i, 0)], tmp_part_size * y_points, MPI_DOUBLE, proc_idx, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                std::string filename = std::string("results/temperature/t=") + std::to_string(current_time) + std::string(".txt");
                write_to_file(filename, temperature, x, y);

                std::string exact_filename = std::string("results/exact_temperature/t=") + std::to_string(current_time) + std::string(".txt");
                write_to_file(exact_filename, exact_temperature, x, y); 
            }                     
        }

        struct timeval write_tv;
        gettimeofday(&write_tv, NULL);

        solution_time += (sol_tv.tv_sec - ref_tv.tv_sec) + (sol_tv.tv_usec - ref_tv.tv_usec)*1.0e-6;
        write_time += (write_tv.tv_sec - sol_tv.tv_sec) + (write_tv.tv_usec - sol_tv.tv_usec)*1.0e-6;
    }

    if (cur_proc == 0) {    
        printf("\nCalculation is running at the %d thread(s)\n", proc_num);
    }

    if (cur_proc == 0) {
        printf("\n[Init phase] execution time = %.6f sec\n", init_time);
        printf("\n[Solution phase] execution time = %.6f sec\n", solution_time);
        printf("\n[Write results phase] execution time = %.6f sec\n", write_time);
    }

    delete[] temperature;

    delete[] exact_temperature;

    MPI_Finalize();

    return 0;
}
