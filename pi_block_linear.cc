#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: init MPI
    srand(time(NULL * world_rank));

    long long int* local_count =(long long int*)malloc(sizeof(long long int) * world_rank);
    long long int total_count;
    if (world_rank > 0)
    {
        // TODO: handle workers
        int number_in_circle = 0;

        for (int i = 0; i < toss / world_rank; i++) {
            float x = fRand(-1, 1);
            float y = fRand(-1, 1);
            float distance_squared = x * x + y * y;
            if (distance_squared <= 1) {
                number_in_circle++;
            }
        }
        MPI_Send(
        /* data         = */ &count, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_LONG_LONG, 
        /* destination  = */ 0, 
        /* tag          = */ 0, 
        /* communicator = */ MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // TODO: master
        for (int i = 0; i < toss / world_rank; i++) {
            float x = fRand(-1, 1);
            float y = fRand(-1, 1);
            float distance_squared = x * x + y * y;
            if (distance_squared <= 1) {
                local_count[0]++;
            }
        }
        for (int i = 1; i< world_rank; i++) {
            MPI_Recv(
            /* data         = */ &local_count[i], 
            /* count        = */ 1, 
            /* datatype     = */ MPI_LONG_LONG, 
            /* source       = */ MPI_ANYSOURCE, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD, 
            /* status       = */ MPI_STATUS_IGNORE);
        }
    }

    if (world_rank == 0)
    {
        // TODO: process PI result
        for (int i = 0; i < world_rank; i++) {
            total_count += local_count[i];
        }
        pi_result = ((double)total_count / (double)toss) * 4.0;
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
