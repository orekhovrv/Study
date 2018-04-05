#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define LENGTH 1.0
#define U_INITIAL 1.0
#define U_BOUNDARY 0.0
#define K_PARAMETER 1.0
#define PARAMETER 0.5
#define POINTS_OUTPUT 11
#define N_EXACT 20



double exact_solution(double x, double t)
{
    unsigned int i;
    double u = 0;
    double x_new, t_new;

    for (i = 0; i <= N_EXACT; ++i)
    {
        x_new = M_PI * (2 * i + 1) * x / LENGTH;
        t_new = (-K_PARAMETER) * (M_PI * M_PI) * (2 * i + 1) * (2 * i + 1) * 
                t / (LENGTH * LENGTH);

        u += exp(t_new) * sin(x_new) / (2 * i + 1);
    }

    u *= 4.0 * U_INITIAL / M_PI;

    return u;
}



int main (int argc, char * argv[])
{
    unsigned int points_count = atoi(argv[1]);
    double final_time = atof(argv[2]);

    unsigned int processes_count, my_rank;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    unsigned int my_points_count;

    unsigned int points_count_0 = points_count / processes_count;
    unsigned int points_count_1 = points_count_0 + 1;

    unsigned int processes_count_1 = 
            points_count - points_count_0 * processes_count;
    unsigned int processes_count_0 = processes_count - processes_count_1;

    if (my_rank < processes_count_0)
        my_points_count = points_count_0;
    else
        my_points_count = points_count_1;

    my_points_count += 2;

    double *u_value = malloc(my_points_count * sizeof(double));

    unsigned int point_number;

    for (point_number = 0; point_number < my_points_count; ++point_number)
        u_value[point_number] = U_INITIAL;

    if (my_rank == 0)
        u_value[0] = U_BOUNDARY;

    if (my_rank == (processes_count - 1))
        u_value[my_points_count - 1] = U_BOUNDARY;

    double x_step = LENGTH / (points_count + 1);
    double t_step = (x_step * x_step * PARAMETER) / K_PARAMETER;

    double current_time;
    double previous_value, current_value;

    double start_time = MPI_Wtime();

    for (current_time = 0.0; current_time < final_time; current_time += t_step)
    {
        previous_value = u_value[0];
        current_value = u_value[1];

        for (point_number = 1; point_number <= my_points_count - 2; 
                ++point_number)
        {
            u_value[point_number] += PARAMETER * (previous_value 
                    - 2 * current_value + u_value[point_number + 1]);

            previous_value = current_value;
            current_value = u_value[point_number + 1];
        }

        if (my_rank % 2 == 0)
        {
            if (my_rank != processes_count - 1)
            {
                MPI_Send(&u_value[my_points_count - 2], 1, MPI_DOUBLE, 
                        my_rank + 1, my_rank + 1, MPI_COMM_WORLD);

                MPI_Recv(&u_value[my_points_count - 1], 1, MPI_DOUBLE, 
                        my_rank + 1, my_rank, MPI_COMM_WORLD, &status);
            }

            if (my_rank != 0)
            {
                MPI_Send(&u_value[1], 1, MPI_DOUBLE, my_rank - 1, my_rank - 1, 
                        MPI_COMM_WORLD);

                MPI_Recv(&u_value[0], 1, MPI_DOUBLE, my_rank - 1, my_rank, 
                        MPI_COMM_WORLD, &status);
            }
        }
        else
        {
            MPI_Recv(&u_value[0], 1, MPI_DOUBLE, my_rank - 1, my_rank, 
                    MPI_COMM_WORLD, &status);

            MPI_Send(&u_value[1], 1, MPI_DOUBLE, my_rank - 1, my_rank - 1, 
                    MPI_COMM_WORLD);

            if (my_rank != processes_count - 1)
            {
                MPI_Recv(&u_value[my_points_count - 1], 1, MPI_DOUBLE, 
                        my_rank + 1, my_rank, MPI_COMM_WORLD, &status);

                MPI_Send(&u_value[my_points_count - 2], 1, MPI_DOUBLE, 
                        my_rank + 1, my_rank + 1, MPI_COMM_WORLD);
            }
        }
    }

    unsigned int previous_points;

    if (my_rank < processes_count_0)
        previous_points = my_rank * points_count_0;
    else
        previous_points = processes_count_0 * points_count_0 + 
                (my_rank - processes_count_0) * points_count_1;

    unsigned int first_point = previous_points + 1;
    unsigned int last_point = previous_points + my_points_count - 2;

    unsigned int points_output_step_0 = points_count / (POINTS_OUTPUT - 1);
    unsigned int points_output_step_1 = points_output_step_0 + 1;

    unsigned int points_output_step_1_steps = 
            points_count - (POINTS_OUTPUT - 1) * points_output_step_0;
    unsigned int points_output_step_0_steps = 
            POINTS_OUTPUT - 2 - points_output_step_1_steps;

    char continue_flag = '0';

    if (my_rank == 0)
    {
        printf("processes count = %d points count = %d final time = %f \n", 
                processes_count, points_count, final_time);
//        printf("u[%f] = %f \n", 0.0, u_value[0]);

        continue_flag = '1';
    }
    else
        MPI_Recv(&continue_flag, 1, MPI_CHAR, my_rank - 1, my_rank, 
                MPI_COMM_WORLD, &status);

    unsigned int point_number_output = 0;
    double x;

    for (point_number = 1; point_number < POINTS_OUTPUT - 1; ++point_number)
    {
        if (point_number <= points_output_step_0_steps)
            point_number_output += points_output_step_0;
        else
            point_number_output += points_output_step_1;

        if (first_point <= point_number_output && 
                point_number_output <= last_point)
        {
            x = point_number_output * x_step;

//            printf("u[%f] = %f \n", 
//                    x, u_value[point_number_output - first_point + 1]);
        }
    }

    if (my_rank == processes_count - 1)
    {
//        printf("u[%f] = %f \n", LENGTH, u_value[my_points_count - 1]);

        double execution_time = MPI_Wtime() - start_time;

        printf("execution time = %f \n", execution_time);

        for (point_number = 0; point_number < POINTS_OUTPUT; ++point_number)
        {
            x = LENGTH * point_number / (POINTS_OUTPUT - 1);

//            printf("exact u[%f] = %f \n", x, exact_solution(x, final_time));
        }
    }
    else
        MPI_Send(&continue_flag, 1, MPI_CHAR, my_rank + 1, my_rank + 1, 
                MPI_COMM_WORLD);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
