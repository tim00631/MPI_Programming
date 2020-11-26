#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <stdint.h>

struct xorshift128p_state {
  uint64_t a, b;
};

/* The state must be seeded so that it is not all zero */
uint64_t xorshift128p(struct xorshift128p_state *state)
{
	uint64_t t = state->a;
	uint64_t const s = state->b;
	state->a = s;
	t ^= t << 23;		// a
	t ^= t >> 17;		// b
	t ^= s ^ (s >> 26);	// c
	state->b = t;
	return t + s;
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

    // TODO: MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    uint32_t seed = time(NULL) * (world_rank + 1);
    struct xorshift128p_state* state = (struct xorshift128p_state*)malloc(sizeof(struct xorshift128p_state));
    state->a = seed + 1;
    state->b = seed & 0x55555555;

    // TODO: binary tree reduction
    
    // ===== pi Estimation Block start =====
    uint64_t number_in_circle = 0;
    uint64_t max_iter = tosses / world_size;
    for (uint64_t i = 0; i < max_iter; i++) {
        uint64_t tmp = xorshift128p(state);
        double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
        double y = (double)(tmp >> 32) / __UINT32_MAX__;
        double distance_squared = x * x + y * y;
        if (distance_squared <= 1) {
            number_in_circle++;
        } 
    }
    // ===== pi Estimation Block end =====

    if (world_rank % 2 == 1) {
        MPI_Send(
        /* data         = */ &number_in_circle, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_LONG_LONG, 
        /* destination  = */ world_rank - 1, 
        /* tag          = */ 0, 
        /* communicator = */ MPI_COMM_WORLD);
    }
    else {
        int s = 1;
        for (int i = 0; i < log(world_size); i++) {
            if(world_rank == 0) {
                uint64_t rcv_temp = 0;
                MPI_Recv(
                /* data         = */ &rcv_temp, 
                /* count        = */ 1, 
                /* datatype     = */ MPI_LONG_LONG, 
                /* source       = */ world_rank + s, 
                /* tag          = */ 0,
                /* communicator = */ MPI_COMM_WORLD, 
                /* status       = */ MPI_STATUS_IGNORE);  
                number_in_circle += rcv_temp; 
            }
            else {
                if (world_rank + s < world_size && world_rank - s > 0) {
                    uint64_t rcv_temp = 0;
                    MPI_Recv(
                    /* data         = */ &rcv_temp, 
                    /* count        = */ 1, 
                    /* datatype     = */ MPI_LONG_LONG, 
                    /* source       = */ world_rank + s, 
                    /* tag          = */ 0,
                    /* communicator = */ MPI_COMM_WORLD, 
                    /* status       = */ MPI_STATUS_IGNORE);  
                    number_in_circle += rcv_temp;
                }
                else {
                    MPI_Send(
                    /* data         = */ number_in_circle, 
                    /* count        = */ 1, 
                    /* datatype     = */ MPI_LONG_LONG, 
                    /* source       = */ world_rank - s, 
                    /* tag          = */ 0,
                    /* communicator = */ MPI_COMM_WORLD); 
                }
            }
            s = s * 2;
        }
    }

    if (world_rank == 0) {
        // TODO: PI result
        pi_result = ((double)number_in_circle / (double)tosses) * 4.0;
        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
