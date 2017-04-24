#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "common.h"
#include "omp.h"

// maximum number of partitions for each axis
#define MAX_PARTITIONS 3000
#define VIT vector<int>::iterator

using namespace std;

vector<int> bin[MAX_PARTITIONS][MAX_PARTITIONS];
omp_lock_t locks[MAX_PARTITIONS][MAX_PARTITIONS]; 

//
//  main program
//
int main(int argc, char **argv) {   
    int navg, nabsavg=0, numthreads;
    double dmin, absmin=1.0, davg, absavg=0.0;

    if(find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n"); 
        printf("-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen (sumname, "a") : NULL;      

    set_size(n);
    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
    init_particles(n, particles);

    const int num_partitions = MAX_PARTITIONS;
    const double size = get_size();
    const double partition_size = max(cutoff, (double) size / num_partitions);
    for (int i = 0 ; i < n ; i ++) {
        int px = (int)(particles[i].x / partition_size);
        int py = (int)(particles[i].y / partition_size);
        bin[px][py].push_back(i);
    }
    
    for (int i = 0 ; i < num_partitions ; i++) {
 	  for (int j = 0 ; j < num_partitions; j++) omp_init_lock(&(locks[i][j]));
    }

    double simulation_time = read_timer();

    #pragma omp parallel private(dmin) 
    {
         numthreads = omp_get_num_threads(); 
    for (int step = 0 ; step < NSTEPS ; step ++) {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        
        //
        //  compute all forces
        //
        #pragma omp for reduction(+:navg) reduction(+:davg)
        for(int i = 0; i < n ; i ++) {
            particles[i].ax = particles[i].ay = 0;
            int px = (int)(particles[i].x / partition_size);
            int py = (int)(particles[i].y / partition_size);
            // only need to compute the forces exerted by particle inside the same cell
            // and the surrounding cells.
            for (int j = max(0, px - 1) ; j <= min(px + 1, num_partitions - 1) ; j ++) {
                for (int k = max(0, py - 1) ; k <= min(py + 1, num_partitions - 1) ; k ++) {
                    if (j == px && k == py && bin[j][k].size() == 1) continue;
                    for (VIT it = bin[j][k].begin() ; it != bin[j][k].end() ; it ++) {
                        apply_force(particles[i], particles[*it], &dmin, &davg, &navg);
                    }
                }
            }
        }
  
        // move particles and update each particle's partition
	    #pragma omp for
        for (int i = 0 ; i < n ; i ++) {
            int ix = (int)(particles[i].x / partition_size);
            int iy = (int)(particles[i].y / partition_size);
            move(particles[i]);
            int px = (int)(particles[i].x / partition_size);
            int py = (int)(particles[i].y / partition_size);
             
            if (ix != px || iy != py) {
                // remove i from the old cell
                omp_set_lock(&(locks[ix][iy]));
                for (VIT it = bin[ix][iy].begin() ; it != bin[ix][iy].end() ; it ++) {
                    if (*it == i) {
                        *it = bin[ix][iy].back();
                        break;
                    }
                }
                bin[ix][iy].pop_back();
                omp_unset_lock(&(locks[ix][iy]));

                // add i to the new cell
                omp_set_lock(&(locks[px][py]));
                bin[px][py].push_back(i);
		        omp_unset_lock(&(locks[px][py]));
            }
        }

        if (find_option( argc, argv, "-no" ) == -1) {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
          if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if (fsave && (step % SAVEFREQ) == 0) save(fsave, n, particles);
        }
    }
    }
    simulation_time = read_timer() - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n, numthreads, simulation_time);

    if (find_option( argc, argv, "-no" ) == -1) {
        if (nabsavg) absavg /= nabsavg;
        // 
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf(", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if (fsum) fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if (fsum) fclose(fsum);
    free(particles);
    if (fsave) fclose(fsave);
    
    return 0;
}
