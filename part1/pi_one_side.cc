#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>

// struct xorshift128p_state {
//   uint64_t a, b;
// };

// /* The state must be seeded so that it is not all zero */
// uint64_t xorshift128p(struct xorshift128p_state *state)
// {
// 	uint64_t t = state->a;
// 	uint64_t const s = state->b;
// 	state->a = s;
// 	t ^= t << 23;		// a
// 	t ^= t >> 17;		// b
// 	t ^= s ^ (s >> 26);	// c
// 	state->b = t;
// 	return t + s;
// }
uint64_t xorshift128p(uint64_t *s)
{
	uint64_t t = s[0];
	uint64_t const r = s[1];
	s[0] = r;
	t ^= t << 23;		// a
	t ^= t >> 17;		// b
	t ^= r ^ (r >> 26);	// c
	s[1]= t;
	return t + r;
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

    MPI_Win win;
    

    // TODO: MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    uint32_t seed = time(NULL) * (world_rank + 1);
    uint64_t s[2];
    s[0] = seed;
    s[1] = seed & 0x55555555;
    // struct xorshift128p_state* state = (struct xorshift128p_state*)malloc(sizeof(struct xorshift128p_state));
    // state->a = seed + 1;
    // state->b = seed & 0x55555555;
    uint64_t* local_count;
    // ===== pi Estimation Block start =====
    uint64_t number_in_circle = 0;
    uint64_t max_iter = tosses / world_size;
    for (uint64_t i = 0; i < max_iter; i++) {
        uint64_t tmp = xorshift128p(s);
        double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
        double y = (double)(tmp >> 32) / __UINT32_MAX__;
        double distance_squared = x * x + y * y;
        if (distance_squared <= 1) {
            number_in_circle++;
        } 
    }
    if (world_rank == 0)
    {
        // total_count += number_in_circle;
        // Master
        local_count = (uint64_t *) malloc(sizeof(uint64_t) * world_size);
        MPI_Win_create(local_count, sizeof(uint64_t) * world_size, sizeof(uint64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        MPI_Put(&number_in_circle, 1, MPI_LONG_LONG, 0, world_rank, 1, MPI_LONG_LONG, win);
        MPI_Win_unlock(0, win);
    }
    else
    {
        // Workers
        MPI_Win_create(local_count, 0, sizeof(uint64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        // MPI_Get(&total_count, 1, MPI_LONG_LONG, 0, 0, 1, MPI_LONG_LONG, win);
        // total_count += number_in_circle;
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        MPI_Put(&number_in_circle, 1, MPI_LONG_LONG, 0, world_rank, 1, MPI_LONG_LONG, win);
        MPI_Win_unlock(0, win);
        
    }

    MPI_Win_free(&win);

    if (world_rank == 0)
    {
        // TODO: handle PI result
        uint64_t total_count = 0;

        for (int i = 0; i < world_size; i++) {
            // printf("local_count[%d]:%lld\n", i, local_count[i]);
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