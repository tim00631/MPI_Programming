#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <immintrin.h>

#define UINT_32_MAX 0xFFFFFFFF
// __uint32_t xor128(__uint32_t seed) {
//     static __uint32_t x = seed;
//     static __uint32_t y = 362436069;
//     static __uint32_t z = 521288629;
//     static __uint32_t w = 88675123;
//     __uint32_t t;

//     t = x ^ (x << 11);
//     x = y;
//     y = z;
//     z = w;
//     return w = w ^ (w >> 19) ^ t ^ (t >> 8);
// }



// double fRand(__uint32_t seed) {
//     // long long MAX = ((long long)RAND_MAX << 31) + RAND_MAX;
//     // long long rand_num = ((long long)rand() << 31) + rand();
//     // printf("%lld, %lld, %lf\n", rand_num, MAX, (double)rand_num/MAX);

//     return xor128(seed) / 4294967296.0;
//     // return ((double)rand_num / (double)MAX);
// }

struct xorshift128_state {
  u_int32_t a, b, c, d;
};

/* The state array must be initialized to not be all zero */
u_int32_t xorshift128(struct xorshift128_state *state)
{
	/* Algorithm "xor128" from p. 5 of Marsaglia, "Xorshift RNGs" */
	u_int32_t t = state->d;

	u_int32_t const s = state->a;
	state->d = state->c;
	state->c = state->b;
	state->b = s;

	t ^= t << 11;
	t ^= t >> 8;
	return state->a = t ^ s ^ (s >> 19);
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
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    u_int32_t seed = time(NULL) * world_rank;
    struct xorshift128_state* state = (struct xorshift128_state*)malloc(sizeof(struct xorshift128_state));
    state->a = seed;
    state->b = seed << 1;
    state->c = seed << 2;
    state->d = seed >> 4;
    long long iteration = tosses / world_size; 
    long long int* local_count;

    if (world_rank > 0) {
        // TODO: handle workers
        long long int number_in_circle = 0;

        for (int i = 0; i < iteration; i++) {
            float x = (float)xorshift128(state) / UINT_32_MAX;
            float y = (float)xorshift128(state) / UINT_32_MAX;
            float distance_squared = x * x + y * y;
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
            float x = xorshift128(state)/UINT_32_MAX;
            float y = xorshift128(state)/UINT_32_MAX;
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
