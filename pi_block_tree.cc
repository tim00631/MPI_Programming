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

    uint64_t number_in_circle = 0;

    
    // TODO: binary tree reduction

    // ===== pi Estimation Block start =====
    // for (__uint64_t i = 0; i < tosses / world_size; i++) {
    //     uint64_t tmp = xorshift128p(state);
    //     double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
    //     double y = (double)(tmp >> 32) / __UINT32_MAX__;
    //     float distance_squared = x * x + y * y;
    //     if (distance_squared <= 1) {
    //         number_in_circle++;
    //     }
    // }
    // ===== pi Estimation Block end =====

    // int round = log2(world_size);
    // int diff = 1;
    // int jump = 2;
    // for (int iter = 0; i < round; iter++) 
    // {            
    //     for (int i = 0; i< world_size; i+=diff) 
    //     {

    //     }
    //     if(world_rank % jump == 1)
    //     {
    //         MPI_Send(
    //         /* data         = */ &number_in_circle, 
    //         /* count        = */ 1, 
    //         /* datatype     = */ MPI_LONG_LONG, 
    //         /* destination  = */ world_rank - diff, 
    //         /* tag          = */ 0, 
    //         /* communicator = */ MPI_COMM_WORLD);
    //         break;
    //     }
    //     else {
    //         MPI_Recv(
    //         /* data         = */ &local_count[i], 
    //         /* count        = */ 1, 
    //         /* datatype     = */ MPI_LONG_LONG, 
    //         /* source       = */ world_rank + diff,
    //         /* tag          = */ 0,
    //         /* communicator = */ MPI_COMM_WORLD, 
    //         /* status       = */ MPI_STATUS_IGNORE);
    //     }
    //     diff = jump;
    //     jump = jump << 1;
    // }
    if (world_rank == 0)
    {
        // TODO: PI result

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
