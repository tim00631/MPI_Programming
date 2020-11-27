#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
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
    // TODO: use MPI_Reduce
    uint32_t seed = time(NULL) * (world_rank + 1);
    struct xorshift128p_state* state = (struct xorshift128p_state*)malloc(sizeof(struct xorshift128p_state));
    state->a = seed + 1;
    state->b = seed & 0x55555555;

    uint64_t* total_count;
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
    if (world_rank == 0) {
        total_count = (uint64_t*)malloc(sizeof(uint64_t) * world_size);
    }
    MPI_Reduce(
        /* send_data     = */ &number_in_circle, 
        /* recv_data     = */ total_count, 
        /* count         = */ 1,
        /* datatype      = */ MPI_LONG_LONG, 
        /* op            = */ MPI_SUM,
        /* root          = */ 0,
        /* communicator  = */ MPI_COMM_WORLD);
    // ===== pi Estimation Block end =====
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
