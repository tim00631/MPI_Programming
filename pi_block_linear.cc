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

// static inline uint64_t rotl(const uint64_t x, int k) {
// 	return (x << k) | (x >> (64 - k));
// }


// static uint64_t s[2];

// uint64_t next(void) {
// 	const uint64_t s0 = s[0];
// 	uint64_t s1 = s[1];
// 	const uint64_t result = s0 + s1;

// 	s1 ^= s0;
// 	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
// 	s[1] = rotl(s1, 37); // c

// 	return result;
// }

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
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    __uint32_t seed = time(NULL) * (world_rank + 1);

    struct xorshift128p_state* state = (struct xorshift128p_state*)malloc(sizeof(struct xorshift128p_state));
    state->a = seed + 1;
    state->b = seed & 0x55555555;
    // s[0] = seed;
    // s[1] = s[0] & 0x55555555;
    long long iteration = tosses / world_size; 
    long long int* local_count;

    if (world_rank > 0) {
        // TODO: handle workers
        long long int number_in_circle = 0;
        for (int i = 0; i < iteration; i++) {
            // uint64_t tmp = next();
            uint64_t tmp = xorshift128p(state);
            double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
            double y = (double)(tmp >> 32) / __UINT32_MAX__;
            // double x = (double) s[0] / __UINT64_MAX__;
            // double y = (double) s[1] / __UINT64_MAX__;
            double distance_squared = x * x + y * y;
            if (distance_squared <= 1) {
                number_in_circle++;
            }
        }

        MPI_Send(
        /* data         = */ &number_in_circle, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_LONG_LONG, 
        /* destination  = */ 0, 
        /* tag          = */ 0, 
        /* communicator = */ MPI_COMM_WORLD);
    }

    else if (world_rank == 0) {
        // TODO: master
        local_count =(long long int*)malloc(sizeof(long long int) * world_size); // initialize global variable
        long long int number_in_circle = 0;
        for (int i = 0; i < tosses / world_size; i++) {
            // uint64_t tmp = next();
            uint64_t tmp = xorshift128p(state);
            double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
            double y = (double)(tmp >> 32) / __UINT32_MAX__;
            // double x = (double) s[0] / __UINT64_MAX__;
            // double y = (double) s[1] / __UINT64_MAX__;
            float distance_squared = x * x + y * y;
            if (distance_squared <= 1) {
                number_in_circle++;
            }
        }
        local_count[0] = number_in_circle;

        for (int i = 1; i< world_size; i++) {
            MPI_Recv(
            /* data         = */ &local_count[i], 
            /* count        = */ 1, 
            /* datatype     = */ MPI_LONG_LONG, 
            /* source       = */ i, 
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD, 
            /* status       = */ MPI_STATUS_IGNORE);
        }
    }

    if (world_rank == 0) {
        // TODO: process PI result
        long long int total_count = 0;

        for (int i = 0; i < world_size; i++) {
            printf("local_count[%d]:%lld\n", i, local_count[i]);
            total_count += local_count[i];
        }
        pi_result = ((double)total_count / (double)tosses) * 4.0;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
