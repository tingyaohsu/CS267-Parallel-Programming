   __m128d  c0010, c2030, c0111, c2131, c0212, c2232, c0313, c2333;

  register double* c0111_ptr = c + lda;
  register double* c0212_ptr = c0111_ptr + lda;
  register double* c0313_ptr = c0212_ptr + lda;

  c0010 = _mm_loadu_pd(c);
  c2030 = _mm_loadu_pd(c+2);
  c0111 = _mm_loadu_pd(c0111_ptr);
  c2131 = _mm_loadu_pd(c0111_ptr+2);
  c0212 = _mm_loadu_pd(c0212_ptr);
  c2232 = _mm_loadu_pd(c0212_ptr+2);
  c0313 = _mm_loadu_pd(c0313_ptr);
  c2333 = _mm_loadu_pd(c0313_ptr+2);


  __m128d a01, a23, b0, b1, b2, b3;
  for (int i = 0; i < K; i++) 
  {
    a01 = _mm_load_pd(a);
    a23 = _mm_load_pd(a+2);
    a += 4;

    b0 = _mm_loaddup_pd(b++);
    b1 = _mm_loaddup_pd(b++);
    b2 = _mm_loaddup_pd(b++);
    b3 = _mm_loaddup_pd(b++);

    c0010 = _mm_add_pd(c0010, _mm_mul_pd(a01, b0));
    c2030 = _mm_add_pd(c2030, _mm_mul_pd(a23, b0));
    c0111 = _mm_add_pd(c0111, _mm_mul_pd(a01, b1));
    c2131 = _mm_add_pd(c2131, _mm_mul_pd(a23, b1));
    c0212 = _mm_add_pd(c0212, _mm_mul_pd(a01, b2));
    c2232 = _mm_add_pd(c2232, _mm_mul_pd(a23, b2));
    c0313 = _mm_add_pd(c0313, _mm_mul_pd(a01, b3));
    c2333 = _mm_add_pd(c2333, _mm_mul_pd(a23, b3));
  }

  _mm_storeu_pd(c, c0010);
  _mm_storeu_pd((c+2), c2030);
  _mm_storeu_pd(c0111_ptr, c0111);
  _mm_storeu_pd((c0111_ptr+2), c2131);
  _mm_storeu_pd(c0212_ptr, c0212);
  _mm_storeu_pd((c0212_ptr+2), c2232);
  _mm_storeu_pd(c0313_ptr, c0313);
  _mm_storeu_pd((c0313_ptr+2), c2333);