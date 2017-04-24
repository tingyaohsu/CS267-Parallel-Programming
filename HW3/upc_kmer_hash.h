#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>
#include <upc_strict.h>
#include "contig_generation.h"

/* K-mer data structure for UPC*/
typedef struct kmer_t_upc kmer_t_upc;
struct kmer_t_upc{
   char kmer[KMER_PACKED_LENGTH];
   char l_ext;
   char r_ext;
};


/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
char lookup_kmer_upc(shared kmer_t_upc* memory_heap, shared int64_t* next_index, 
                       shared int64_t* hash_table, int64_t hashtable_size, 
                       const unsigned char *kmer)
{
   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable_size, (char*) packedKmer);
   int64_t pos = hash_table[hashval];
   kmer_t_upc tmp_kmer; 
   for (; pos != -1; ) {
/*      upc_memget(dst, src, size)	copy from shared memory to private memory */
      upc_memget(&tmp_kmer, &memory_heap[pos], sizeof(kmer_t_upc));
      if ( memcmp(packedKmer, tmp_kmer.kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
         return tmp_kmer.r_ext;
      }
      pos = next_index[pos];
   }
   
   return 0;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer_upc(shared int64_t* hash_table,  
             shared kmer_t_upc * memory_heap, 
             const unsigned char *kmer, char left_ext, char right_ext, 
             shared int64_t* next_index, int64_t global_kmer, int hashtable_size)
{
   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable_size, (char*) packedKmer);
   kmer_t_upc tmp_kmer;
/*   upc_memcpy(dst, src, size)     copy from shared memory to shared memory */
   memcpy(tmp_kmer.kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
   tmp_kmer.l_ext = left_ext;
   tmp_kmer.r_ext = right_ext;
/*  upc_memput(dst, src, size)    copy from private memory to shared memory */
   upc_memput(&memory_heap[global_kmer], &tmp_kmer, sizeof(kmer_t_upc));
   
   /* Storing global_kmer to the end of of the list at hashtable[hashval] */
   int64_t ptr;
   /* Atomic operations with a pointer-to-shared */
   ptr = bupc_atomicI64_cswap_strict((shared void*)(&hash_table[hashval]), -1, global_kmer);
   while(ptr != -1) {
     ptr = bupc_atomicI64_cswap_strict((shared void*)(&next_index[ptr]), -1, global_kmer);
   }
   
/*   Deallocation functions 
int dealloc_heap(memory_heap_t *memory_heap)
{
   free(memory_heap->heap);
   return 0;
}

int dealloc_hashtable(hash_table_t *hashtable)
{
   free(hashtable->table);
   return 0;
}
*/   
   
   
   return 0;
}

#endif // KMER_HASH_UPC_H
