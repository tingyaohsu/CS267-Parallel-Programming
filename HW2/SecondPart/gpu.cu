#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define NUM_THREADS 256
#define bin_size (cutoff)
// We can prove, assuming the correctness checker (based on minimum distance) is correct,
// that the number of particles in each bin at any given time step is at most 6.
// Please see our report for more details.
#define MAX_PARTICLES_PER_BIN 6

extern double size;

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;

}

__global__ void move_gpu (particle_t * particles, int n, double size)
{

    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;

    particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}

// Assign each particle to correct bin
__global__ void binning_gpu (particle_t *particles, int n, particle_t **bins, 
                             int *num_particles_per_bins, int nbins_side) {
   // Get thread ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= n) return;

    // Assigning bins to threads
    int bin_x_id = floor(particles[tid].x / bin_size);
    int bin_y_id = floor(particles[tid].y / bin_size);
    int bin_id = nbins_side * bin_x_id  + bin_y_id;

    // Indexing without race condition
    int thread_index_in_bin = atomicAdd(num_particles_per_bins+bin_id, 1);
    bins[MAX_PARTICLES_PER_BIN * bin_id + thread_index_in_bin] = particles + tid;
}

__global__ void compute_forces_gpu(particle_t *particles, particle_t **bins,
                                   int *num_particles_per_bins, int nbins_side, int n) {
  // Get particle ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= n) return;

  int bin_x_id = floor(particles[tid].x / bin_size);
  int bin_y_id = floor(particles[tid].y / bin_size);

  // Compute forces
  particles[tid].ax = particles[tid].ay = 0;
  for (int i = max(bin_x_id - 1, 0) ; i < min(bin_x_id + 2, nbins_side) ; i ++) {
        for (int j = max(bin_y_id - 1, 0) ; j < min(bin_y_id + 2, nbins_side) ; j ++) {
            int bid = nbins_side * i + j;
            for (int k = 0 ; k < num_particles_per_bins[bid] ; k ++) {
                apply_force_gpu(particles[tid], *(bins[bid * MAX_PARTICLES_PER_BIN + k]));
            }
        }
    }
}

int main(int argc, char **argv) {    
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        return 0;
    }
    
    int n = read_int(argc, argv, "-n", 1000);

    char *savename = read_string(argc, argv, "-o", NULL);
    
    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    int nbins_side = ceil(size / bin_size);
    int numbins = nbins_side * nbins_side;
    particle_t **d_bins; 
    cudaMalloc((void ***) &d_bins, numbins * MAX_PARTICLES_PER_BIN * sizeof(particle_t *));

    int *num_particles_per_bins;
    cudaMalloc((void **) &num_particles_per_bins, numbins * sizeof(int *));

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer() - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer();

    for( int step = 0; step < NSTEPS; step++ ) {
        //
        //  compute forces
        //
	    int blks = (n + NUM_THREADS - 1) / NUM_THREADS;

        cudaMemset(num_particles_per_bins, 0, numbins * sizeof(int));

		// assign particles to bins
        binning_gpu<<<blks, NUM_THREADS>>>(d_particles, n, d_bins, num_particles_per_bins, nbins_side);

        // compute the forces
        compute_forces_gpu<<<blks, NUM_THREADS>>>(d_particles, d_bins, num_particles_per_bins, nbins_side, n);
        
        // move particles
        move_gpu<<<blks, NUM_THREADS>>>(d_particles, n, size);
        
        //
        //  save if necessary
        //
        if(fsave && (step%SAVEFREQ) == 0) {
            // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save(fsave, n, particles);
        }
    }
    cudaThreadSynchronize();
    simulation_time = read_timer() - simulation_time;
    
    printf("CPU-GPU copy time = %g seconds\n", copy_time);
    printf("n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    
    free(particles);
    cudaFree(d_particles);
    cudaFree(d_bins);
    cudaFree(num_particles_per_bins);
    if(fsave) fclose(fsave);
    
    return 0;
}
