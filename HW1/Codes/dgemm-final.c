#include <mmintrin.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <emmintrin.h>
#include <string.h>

const char* dgemm_desc = "Simple blocked dgemm.";
// Tunned this numbers is shown in another subroutine
#define BLOCK_SIZE 132
// Tunned these numbers is shown in another subroutine
#define L1_Block 400
#define L2_Block 7500

#define min(a,b) (((a)<(b))?(a):(b))

/* Vectorization using SIMD-SSE - 4-by-4 blocks */

static inline void pack_a (int lda, const int K, double* a_unpack_ptr, double* a_pack_ptr) {
  for (int n = 0; n < K; ++n)
  {
    *a_pack_ptr++ = *a_unpack_ptr;
    *a_pack_ptr++ = *(a_unpack_ptr+1);
    *a_pack_ptr++ = *(a_unpack_ptr+2);
    *a_pack_ptr++ = *(a_unpack_ptr+3);
    a_unpack_ptr += lda;
  }
}

static inline void pack_b (int lda, const int K, double*  b_unpack_ptr, double* b_pack_ptr) {
  double *b_unpack_ptr0, *b_unpack_ptr1, *b_unpack_ptr2, *b_unpack_ptr3;
  b_unpack_ptr0 = b_unpack_ptr;
  b_unpack_ptr1 = b_unpack_ptr0 + lda;
  b_unpack_ptr2 = b_unpack_ptr1 + lda;
  b_unpack_ptr3 = b_unpack_ptr2 + lda;

  for (int i = 0; i < K; ++i)
  {
    *b_pack_ptr++ = *b_unpack_ptr0++;
    *b_pack_ptr++ = *b_unpack_ptr1++;
    *b_pack_ptr++ = *b_unpack_ptr2++;
    *b_pack_ptr++ = *b_unpack_ptr3++;
  }
}

static inline void SSE(int lda, int K, double* restrict a, double* restrict b, double* restrict c)
{
#include "SSE.h"
}


/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static inline void do_block (int lda, int M, int N, int K, double* restrict A, double* restrict B, double* restrict C)
{

    double *a_ptr, *b_ptr, *c;
    double Ab[M*K], Bb[K*N];
    int i = 0, j = 0, p = 0;
    for (j = 0 ; j < N-3; j += 4) 
    {
      b_ptr = &Bb[j*K];
      // copy optimization and Bb[] gets transposed
      pack_b(lda, K, B+j*lda, b_ptr);
      for (i = 0; i < M-3; i += 4) {
        a_ptr = &Ab[i*K];
        if (j == 0) pack_a(lda, K, A+i, a_ptr);
        c = C + i + j*lda;
        SSE(lda, K, a_ptr, b_ptr, c);
      }
    }

	// ------------Padding -----------------------
	// Pad if M and N are not multiples of 4
    int mod1 = M%4;
    int mod2 = N%4;
  if (mod1 != 0) 
  {
    for ( ; i < M; ++i)
      for (p = 0; p < N; ++p) 
      {
        /* Compute C[i,j] */
        double c_ip = C[i+p*lda];
        for (int k = 0; k < K; ++k)
          c_ip += A[i+k*lda] * B[k+p*lda];
          C[i+p*lda] = c_ip;
      }
  }
  if (mod2 != 0) 
  {
    for ( ; j < N; ++j)
      for (i = 0; i < (M - mod1); ++i) 
      {
        /* Compute C[i,j] */
        double c_ij = C[i+j*lda];
        for (int k = 0; k < K; ++k)
          c_ij += A[i+k*lda] * B[k+j*lda];
          C[i+j*lda] = c_ij;
      }
  }
	// ------------------------------
}  

// /* This routine performs a dgemm operation
//  *  C := C + A * B
//  * where A, B, and C are lda-by-lda matrices stored in column-major format.
//  * On exit, A and B maintain their input values.
//  * Optimization: Two levels of blocking. */

void square_dgemm (int lda, double* A, double* B, double* C)
{	
	  // Multi level Cache Blocking for L1 and L2
	 // L1_Block and L2_Block are tuned separately and best values are used here
	int l,m,n,i,j,k,ii,jj,kk,K,M,N;
	for (l = 0; l < lda; l += L2_Block) 
	{
		kk = l + min(L2_Block, lda-l);
		/* For each block-column of B of size L2*/
		for (m = 0; m < lda; m += L2_Block) 
		{
			jj = m + min(L2_Block, lda-m);
			/* For each block-row of A of size L2*/ 
			for (n = 0; n < lda; n += L2_Block) 
			{
				ii = n + min(L2_Block, lda-n);
				for (k = l; k < kk; k += L1_Block) 
				{
					K = min(L1_Block, kk-k);
					/* For each block-column of B of size L1*/
					for (j = m; j < jj; j += L1_Block) 
					{
						N = min(L1_Block, jj-j);
						/* For each block-row of A of size L1*/ 
						for (i = n; i < ii; i += L1_Block) 
						{
							M = min(L1_Block, ii-i);
							/* Performs a smaller dgemm operation
							*  C' := C' + A' * B'
							* where C' is M-by-N, A' is M-by-K, and B' is K-by-N. */
							do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
						}
					}
				}
			}
		}
	}
}

