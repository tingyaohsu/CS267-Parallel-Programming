#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>
#include "packingDNAseq.h"
#include "upc_kmer_hash.h"


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	// -----------------------------------------------------
    /* Read the input file name */
    char *input_UFX_name;
    input_UFX_name = argv[1];

    /* Extract the number of k-mers in the input file */
	int64_t nKmers;
        nKmers = getNumKmersInUFX(input_UFX_name);

    /* Split the work of reading kmers between threads and store them in local working_buffer */
	int64_t lines_to_read_per_thread = nKmers / THREADS + (MYTHREAD < (nKmers % THREADS));
	int64_t total_chars_to_read_per_thread = lines_to_read_per_thread * LINE_SIZE;
	int64_t remaining_lines = MYTHREAD * nKmers / THREADS;
	if (MYTHREAD <= nKmers % THREADS) remaining_lines += MYTHREAD;
	else  remaining_lines += nKmers % THREADS;
	int64_t total_remaining_chars = remaining_lines * LINE_SIZE;
	int64_t cur_chars_read;
	upc_file_t *fid;
	unsigned char* private_buffer = (unsigned char*) malloc((total_chars_to_read_per_thread+5) * sizeof(unsigned char));
	
	/* upc_all_fopen : a collective read operation using individual file pointers. Each thread reads a block of data into 
					   a private buffer from a particular thread-specific offset. */
	fid = upc_all_fopen(input_UFX_name, UPC_RDONLY | UPC_INDIVIDUAL_FP, 0, NULL ); 
	
	/* upc_all_fseek : sets the current position of the file pointer associated with fid. */
	upc_all_fseek( fid, total_remaining_chars*sizeof(unsigned char), UPC_SEEK_SET); 
	
	/* upc_all_fread_local : reads data from a file into a local buffer on each thread. */
	cur_chars_read = upc_all_fread_local(fid, private_buffer, sizeof(unsigned char), total_chars_to_read_per_thread,UPC_IN_ALLSYNC | 					UPC_OUT_ALLSYNC);
	upc_barrier;
	upc_all_fclose(fid);
	
	// -----------------------------------------------------
	
	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
    /* Initialize lookup table that will be used for the DNA packing routines */
    init_LookupTable();	
	
 	/* Create heap memory, hash table and initialization */
	/* ---------------- Data structures used here: -----------------------
	   upc_all_alloc(n,b) : collectively allocates nxb bytes of shared data distributed across the threads with a block size of b bytes. It is    							intended to be called by all the threads. It returns a single ptr to the shared allocated space
	   upc_forall : iteration distribution across the threads will take place in round-robin fashion */
	
    shared int64_t* next_index = (shared int64_t*) upc_all_alloc(nKmers, sizeof(int64_t));
    shared kmer_t_upc* memory_heap = (shared kmer_t_upc*) upc_all_alloc(nKmers, sizeof(kmer_t_upc)); 
    int64_t hashtable_size = nKmers * LOAD_FACTOR;
    shared int64_t* hash_table = (shared int64_t*) upc_all_alloc(hashtable_size, sizeof(int64_t));
    int64_t i, ptr=0, global_kmer = remaining_lines; 
   // initial hash table
   upc_forall(i=0; i<hashtable_size; i++; &hash_table[i]) 
     hash_table[i] = -1;
   upc_forall(i=0; i<nKmers; i++; &next_index[i]) 
     next_index[i] = -1;
   //syncronization needed here
   upc_barrier; 

 	while (ptr < cur_chars_read) {
     char left_ext = (char) private_buffer[ptr+KMER_LENGTH+1];
     char right_ext = (char) private_buffer[ptr+KMER_LENGTH+2];
	 
     /* Add k-mer to hash table */
     add_kmer_upc(hash_table, memory_heap, &private_buffer[ptr], left_ext, right_ext, next_index, global_kmer, hashtable_size);
     /* Move to the next k-mer in the input buffer */
     ptr += LINE_SIZE;
	 global_kmer++;
 }
    //syncronization needed here
	upc_barrier;
	constrTime += gettime();
	// -----------------------------------------------------
	
	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	// data structures in graph traversal:
	char cur_contig[MAXIMUM_CONTIG_SIZE+1], unpackedKmer[KMER_LENGTH+1]; 
	char right_ext, left_ext;
	char pgen_output[100];
	sprintf(pgen_output, "pgen%d.out", MYTHREAD);
	FILE *upc_output = fopen(pgen_output, "w"); 
	
        i = 0, ptr=0; 
	for (; i < lines_to_read_per_thread; i ++, ptr += LINE_SIZE) {
	  char left_ext = (char) private_buffer[ptr+KMER_LENGTH+1];
	  if (left_ext != 'F') continue;

	  memcpy(cur_contig, &private_buffer[ptr], KMER_LENGTH * sizeof(char));
	  
	  int64_t posInContig = KMER_LENGTH;
      right_ext = (char) private_buffer[ptr+KMER_LENGTH+2];
    
      /* Keep adding bases while not finding a terminal node */	
    while(right_ext != 'F' && right_ext != 0) {
      cur_contig[posInContig ++] = right_ext;
      
      // hash table lookup for right extension
      right_ext = lookup_kmer_upc(memory_heap, next_index, hash_table, hashtable_size,
                                   (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
    }
    
      /* Print the contig */
       cur_contig[posInContig] = '\0';
       fprintf(upc_output, "%s\n", cur_contig);
}

	fclose(upc_output);	
	

	if(MYTHREAD == 0) {
	  upc_free(hash_table);
	  upc_free(memory_heap);
	  upc_free(next_index);
	}
	
	upc_barrier;
	traversalTime += gettime();
	// -----------------------------------------------------
	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
