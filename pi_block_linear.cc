#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#define UINT64_C(c) (c ## ULL)
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

static u_int64_t s[2];

static inline u_int64_t rotl(const u_int64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

u_int64_t next(void) {
	const u_int64_t s0 = s[0];
	u_int64_t s1 = s[1];
	const u_int64_t result = s0 + s1;

	s1 ^= s0;
	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotl(s1, 37); // c

	return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jump(void) {
	static const u_int64_t JUMP[] = { 0xdf900294d8f554a5, 0x170865df4b3201fc };

	u_int64_t s0 = 0;
	u_int64_t s1 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			next();
		}

	s[0] = s0;
	s[1] = s1;
}


/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations. */

void long_jump(void) {
	static const u_int64_t LONG_JUMP[] = { 0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1 };

	u_int64_t s0 = 0;
	u_int64_t s1 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			next();
		}

	s[0] = s0;
	s[1] = s1;
}

// float fRand(__uint32_t seed) {
//     // long long MAX = ((long long)RAND_MAX << 31) + RAND_MAX;
//     // long long rand_num = ((long long)rand() << 31) + rand();
//     // printf("%lld, %lld, %lf\n", rand_num, MAX, (double)rand_num/MAX);

//     return (float)xor128(seed) / __UINT32_MAX__;
//     // return ((double)rand_num / (double)MAX);
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

    __uint32_t seed = time(NULL) * world_rank;
    s[0] = seed;
    s[1] = seed<<2;
    long long iteration = tosses / world_size; 
    long long int* local_count;

    if (world_rank > 0) {
        // TODO: handle workers
        long long int number_in_circle = 0;
        for (int i = 0; i < iteration; i++) {
            long_jump();
            double x = (double)s[0] / __UINT64_MAX__;
            double y = (double)s[1] / __UINT64_MAX__;
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
            long_jump();
            double x = (double)s[0] / __UINT64_MAX__;
            double y = (double)s[1] / __UINT64_MAX__;
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
