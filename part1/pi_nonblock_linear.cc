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
    
    uint32_t seed = time(NULL) * (world_rank + 1);
    struct xorshift128p_state* state = (struct xorshift128p_state*)malloc(sizeof(struct xorshift128p_state));
    state->a = seed + 1;
    state->b = seed & 0x55555555;

    uint64_t* local_count;

    if (world_rank > 0)
    {
        // TODO: MPI workers
        MPI_Request request;

        // ===== pi Estimation Block start =====
        uint64_t number_in_circle = 0;
        uint64_t max_iter = tosses / world_size;
        for (uint64_t i = 0; i < max_iter; i++)
        {
            uint64_t tmp = xorshift128p(state);
            double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
            double y = (double)(tmp >> 32) / __UINT32_MAX__;
            double distance_squared = x * x + y * y;
            if (distance_squared <= 1) {
                number_in_circle++;
            }
            
        }
        // ===== pi Estimation Block end =====
        MPI_Isend(
        /* data         = */ &number_in_circle, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_LONG_LONG, 
        /* destination  = */ 0, 
        /* tag          = */ 0, 
        /* communicator = */ MPI_COMM_WORLD,
        /* request = */      &request);
    }
    else if (world_rank == 0)
    {
        // TODO: non-blocking MPI communication.
        // Use MPI_Irecv, MPI_Wait or MPI_Waitall.
        
        MPI_Request* requests = (MPI_Request *)malloc(sizeof(MPI_Request) * world_size);
        MPI_Status* status = (MPI_Status *)malloc(sizeof(MPI_Status) * world_size);
        local_count = (uint64_t *)malloc(sizeof(uint64_t)* world_size);

        for(int i = 1; i < world_size; i++)
        {
            MPI_Irecv(
            /* data         = */ &local_count[i], 
            /* count        = */ 1, 
            /* datatype     = */ MPI_LONG_LONG, 
            /* source       = */ i, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD,
            /* request = */      &requests[i]);
        }
        
        // ===== pi Estimation Block start =====
        uint64_t number_in_circle = 0;
        uint64_t max_iter = tosses / world_size;
        for (uint64_t i = 0; i < max_iter; i++)
        {
            uint64_t tmp = xorshift128p(state);
            double x = (double)(tmp << 32 >> 32) / __UINT32_MAX__;
            double y = (double)(tmp >> 32) / __UINT32_MAX__;
            float distance_squared = x * x + y * y;
            if (distance_squared <= 1) 
            {
                number_in_circle++;
            }
        }
        local_count[0] = number_in_circle;
        // ===== pi Estimation Block end =====

        MPI_Waitall(world_size-1,requests+1,status+1);
    }

    if (world_rank == 0)
    {
        // TODO: PI result
        uint64_t total_count = 0;

        for (int i = 0; i < world_size; i++)
        {
            // printf("local_count[%d]:%lu\n", i, local_count[i]);
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