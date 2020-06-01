#ifndef STDINCLUDE_H
#define STDINCLUDE_H

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <emmintrin.h>
#include <vector>
#include "pthread.h"
#include "ProgressBar.hpp"

#if defined _WIN32 || defined _WIN64
#else
#include <unistd.h>
#endif

using namespace std;

#if defined _WIN32 || defined _WIN64
#define __WINDOWS__ 1
#endif

enum {Abin=0, Cbin=3, Gbin=1, Tbin=2};
enum InputType {DNA, PROTEIN, UNKNOWN};

typedef __m128i v4ui;

union v4ui_vector {
  v4ui v;
  unsigned int d[4];
  char c[16];
  unsigned long long ul[2];
  };

typedef float distType;

//RAPID_DISK CONSTANT: percentage of free space compared to the current memory usage before resizing and rebuilding the sorted matrix and cache
const int DATA_STRUCTURE_SIZE_TRESSHOLD = 90;
// MEMORY EFFICIENT RAPIDNJ CONSTANT: the minimum percentage of the sorted matrix which must fit in memory.
const distType MIN_SORTED_MATRIX_SIZE = 0.05f;
//Number of buffers used pr. thread for I/O efficient computation of distance matrices
const static int THREAD_ROW_BUFFER_COUNT = 5;

__inline void* malloc_aligned(size_t  size, unsigned int alignment) {
  void *ptr;
#ifdef __WINDOWS__
  ptr = _aligned_malloc(size, alignment);
#elif defined __APPLE__
  //Mac OS X aligns memory to a 16 byte boundary per default
  ptr = malloc(size);
#else
  int ret = posix_memalign(&ptr, alignment, size);
  if(ret != 0){
    ptr = NULL;
}
#endif
  if (ptr == NULL) {
    cerr << "Could not allocate " << size << " bytes of aligned memory" << endl;
    exit(1);
  }
  return ptr;
}

#endif //STDINCLUDE_H
