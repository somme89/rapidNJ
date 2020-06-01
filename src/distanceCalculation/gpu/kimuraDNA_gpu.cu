/*
__device__ void printBitString(unsigned int data){
  for (int j = 0;j < 32 ; j++) {
    if(j % 4 == 0 && j != 0){
      printf(" ");
    }
    printf("%d",((data  >> j) & 1));
  }
  printf("\n");
}


__device__ void printDNASequence(unsigned int data){
  for (int j = 0;j < 16; j++) {
    if(j % 4 == 0 && j != 0){
      printf(" ");
    }
    unsigned int c = ((data >> (j*2)) & 3);
    switch (c) {
    case 0: {
      printf("A ");
      break;
    }
    case 1: {
      printf("G ");
      break;
    }
    case 2: {
      printf("T ");
      break;
    }
    case 3: {
      printf("C ");
      break;
    }	  
    default:
      break;
    }
  }
  printf("\n");
}
*/

#include "cutil_inline.h"
#include "constants.hpp"
#include <sys/time.h>

__shared__ unsigned int shared_data_dna[threads_pr_block*3];

__device__ void parallelReductionDNA2(){
    //wait for all threads to finish and reuse shared memory to store results
  __syncthreads();
  
  //perform parallel reduction
  for(unsigned int s=blockDim.x/2; s>32; s>>=1){    
    if (threadIdx.x < s){
      shared_data_dna[threadIdx.x] += shared_data_dna[threadIdx.x + s];
      shared_data_dna[threadIdx.x + threads_pr_block] += shared_data_dna[threadIdx.x + s + threads_pr_block];
      shared_data_dna[threadIdx.x + threads_pr_block*2] += shared_data_dna[threadIdx.x + s + threads_pr_block*2];
    }
    __syncthreads();
  }
  
  if (threadIdx.x < 32) {
    shared_data_dna[threadIdx.x] += shared_data_dna[threadIdx.x + 32];
    shared_data_dna[threadIdx.x + threads_pr_block] += shared_data_dna[threadIdx.x + 32 + threads_pr_block];
    shared_data_dna[threadIdx.x + threads_pr_block*2] += shared_data_dna[threadIdx.x + 32 + threads_pr_block*2];
  }   
  __syncthreads();
  if (threadIdx.x < 16) {
    shared_data_dna[threadIdx.x] += shared_data_dna[threadIdx.x + 16];
    shared_data_dna[threadIdx.x + threads_pr_block] += shared_data_dna[threadIdx.x + 16 + threads_pr_block];
    shared_data_dna[threadIdx.x + threads_pr_block*2] += shared_data_dna[threadIdx.x + 16 + threads_pr_block*2];
  }
  __syncthreads();
  if (threadIdx.x < 8) {
    shared_data_dna[threadIdx.x] += shared_data_dna[threadIdx.x + 8];
    shared_data_dna[threadIdx.x + threads_pr_block] += shared_data_dna[threadIdx.x +  8 + threads_pr_block];
    shared_data_dna[threadIdx.x + threads_pr_block*2] += shared_data_dna[threadIdx.x +  8 + threads_pr_block*2];
  }
  __syncthreads();
  if (threadIdx.x < 4) {
    shared_data_dna[threadIdx.x] += shared_data_dna[threadIdx.x + 4];
    shared_data_dna[threadIdx.x + threads_pr_block] += shared_data_dna[threadIdx.x +  4 + threads_pr_block];
    shared_data_dna[threadIdx.x + threads_pr_block*2] += shared_data_dna[threadIdx.x +  4 + threads_pr_block*2];
  }  
  __syncthreads();
  
  if (threadIdx.x < 2) {
    shared_data_dna[threadIdx.x] += shared_data_dna[threadIdx.x + 2];
    shared_data_dna[threadIdx.x + threads_pr_block] += shared_data_dna[threadIdx.x + 2 + threads_pr_block];
    shared_data_dna[threadIdx.x + threads_pr_block*2] += shared_data_dna[threadIdx.x + 2 + threads_pr_block*2];
  }  
  __syncthreads();

  if(threadIdx.x == 0){
    shared_data_dna[0] += shared_data_dna[1];
    shared_data_dna[0 + threads_pr_block] += shared_data_dna[1 + threads_pr_block];
    shared_data_dna[0 + threads_pr_block*2] += shared_data_dna[1 + threads_pr_block*2];
  }
  __syncthreads();
}

/*each block computes all distances between a sequence and all other sequences with higher indexes. In this way no summation is needed but shared memory is not utillised*/
__global__ void computeSliceDNA(float* results, unsigned int* bitStrings_gpu, unsigned int* gapFilters_gpu, unsigned int bsStride, unsigned int dataSize, unsigned int rowOffset){

  int i = blockIdx.y;
  if( i < blockIdx.x+1+rowOffset){
    return;
  }
  unsigned int iterations = dataSize / threads_pr_block;
  if(dataSize % threads_pr_block != 0){
    iterations++;
  }
  //for(int i = blockIdx.x+1+rowOffset; i < sequenceCount_GPU; i++) {
  //set counters in shared memory
  shared_data_dna[threadIdx.x] = 0; 
  shared_data_dna[threadIdx.x + threads_pr_block] = 0;
  shared_data_dna[threadIdx.x + threads_pr_block*2] = 0;
  for(int slice = 0; slice < iterations; slice++) {     
    // fetch next block of sequence data
    unsigned int idx0 = (rowOffset + blockIdx.x)*bsStride + slice*threads_pr_block + threadIdx.x; 
    unsigned int baseSeq = bitStrings_gpu[idx0];
    unsigned int baseGap = gapFilters_gpu[idx0];
    unsigned int idx1 = i*bsStride + slice*threads_pr_block + threadIdx.x;
    //printf("%d*%d+%d*%d+%d=%d \n",i,bsStride,slice,threads_pr_block,threadIdx.x,idx1);
    unsigned int r = bitStrings_gpu[idx1];
    unsigned int gf = gapFilters_gpu[idx1];
      
    // compute distances and gaps
    r = r ^ baseSeq;
    gf = gf & baseGap;
    unsigned int tv = (r >> 1) & 0x55555555;
    unsigned int ts = (~tv) & (r & 0x55555555);
      
    // handle gaps
    tv = tv & gf;
    ts = ts & gf;
    // sum distances in this thread
    tv = __popc(tv);
    ts = __popc(ts);
    gf = __popc(gf);
      
    //Remove invalid results
    if(slice * threads_pr_block + threadIdx.x >= dataSize) {
      tv = 0;
      ts = 0;
      gf = 0;
    }
      
    //save intermediate result
    shared_data_dna[threadIdx.x] += tv;
    shared_data_dna[threadIdx.x + threads_pr_block] += ts;
    shared_data_dna[threadIdx.x + threads_pr_block*2] += gf;    
  }

  //sum over all threads
  parallelReductionDNA2();

  //store data in global memory
  if(threadIdx.x == 0){
    //printf("%d \n", blockIdx.x * sequenceCount_GPU + i);
    //tv_gpu[blockIdx.x * gridDim.y + i] = shared_data_dna[0];
    //ts_gpu[blockIdx.x * gridDim.y + i] = shared_data_dna[threads_pr_block];
    unsigned int lengthWithoutGaps = shared_data_dna[threads_pr_block*2];
    float ts = float(shared_data_dna[threads_pr_block]) / lengthWithoutGaps;
    float tv = float(shared_data_dna[0]) / lengthWithoutGaps;
    float temp1 = 1.0f-2.0f*ts-tv;
    float temp2 = 1.0f-2.0f*tv;
    float distance = -1.0f;
    if(!(temp1 <= 0 || temp2 <= 0)){
      distance = -0.5f*log(1.0f-2.0f*ts-tv)-0.25f*log(1.0f-2.0f*tv);
    }
    results[blockIdx.x * gridDim.y + i] = distance;
  }
}

unsigned int sequenceCount2_dna;
unsigned int dataSize2_dna;
size_t resultMemSize2_dna;
size_t bitStringsMemSize2_dna;

//results
float* results_dna_gpu;
unsigned int bsStride2_dna;

//data
unsigned int* gapFilters2_dna_gpu;
unsigned int* bitStrings2_dna_gpu;
unsigned int rowsPrKernel2;
unsigned int numberOfKernelLaunches2;

#ifdef TIMING
extern float totalGpuComputation;
extern timeval start,end;
extern float totalTransfer;
#endif

extern "C" void computeDistancesDNA2_gpu(float* results) {

  // execute the kernel
  //  unsigned int startIdx = i * threads_pr_block;
  for(int i = 0; i < numberOfKernelLaunches2; i++) {
    printf("kernel: %d -----------------------\n",i);
    // setup execution parameters
    unsigned int gridSize_dna = 0;
    if(i != numberOfKernelLaunches2-1){
      gridSize_dna = rowsPrKernel2;
    } else {
      gridSize_dna = sequenceCount2_dna - (numberOfKernelLaunches2-1)*rowsPrKernel2;
    }
    dim3 grid(gridSize_dna,sequenceCount2_dna);
    dim3 block(threads_pr_block, 1, 1);
    unsigned int transferSize = gridSize_dna * sequenceCount2_dna * sizeof(unsigned int);
    computeSliceDNA<<< grid, block >>>(results_dna_gpu, bitStrings2_dna_gpu, gapFilters2_dna_gpu, bsStride2_dna, dataSize2_dna, i*rowsPrKernel2);
    
#ifdef TIMING
    cudaThreadSynchronize();
    gettimeofday(&end,NULL);
    totalGpuComputation += (end.tv_sec - start.tv_sec)*1000.0 + (end.tv_usec - start.tv_usec)/1000.0;
    gettimeofday(&start,NULL);
#endif

    cutilSafeCall(cudaMemcpy(&results[i * rowsPrKernel2 * sequenceCount2_dna], results_dna_gpu, transferSize, cudaMemcpyDeviceToHost));
  
#ifdef TIMING
  cudaThreadSynchronize();
  gettimeofday(&end,NULL);
  totalTransfer += (end.tv_sec - start.tv_sec)*1000.0 + (end.tv_usec - start.tv_usec)/1000.0;
  gettimeofday(&start,NULL);
#endif
  }
}

extern "C" void storeDataDNA2_gpu(unsigned int* bitStrings, unsigned int* gapFilters) {
  // Copy results from device to host
  cutilSafeCall(cudaMemcpyAsync(bitStrings2_dna_gpu, bitStrings, bitStringsMemSize2_dna, cudaMemcpyHostToDevice,0));
  cutilSafeCall(cudaMemcpyAsync(gapFilters2_dna_gpu, gapFilters, bitStringsMemSize2_dna, cudaMemcpyHostToDevice,0));
}

extern "C" void initialiseDNA2_gpu(unsigned int sequenceCount, unsigned int bitStringCount, unsigned int _bsStride) {
  printf("initialising GPU... \n");
  sequenceCount2_dna = sequenceCount;
  bsStride2_dna = _bsStride;
  dataSize2_dna = min(bsStride2_dna,bitStringCount * 4);

  resultMemSize2_dna = sequenceCount * sequenceCount * sizeof(float);
  bitStringsMemSize2_dna = bsStride2_dna * sizeof(unsigned int) * sequenceCount;
  
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);  
  float gpuMemSize = deviceProp.totalGlobalMem * 0.8f;
  
  printf("Gpu memsize %fMB\n",gpuMemSize/1024.0f/1024.0f);

  //printf("TODO REMOVE THIS \n");
  //float gpuMemSize = 50.0 * 1024 * 1024;

  
  rowsPrKernel2 = sequenceCount * min(1.0f,(gpuMemSize - bitStringsMemSize2_dna*2) / resultMemSize2_dna);
  numberOfKernelLaunches2 = sequenceCount/rowsPrKernel2;
  if(sequenceCount % rowsPrKernel2 != 0){
    numberOfKernelLaunches2++;
  }
  rowsPrKernel2 = sequenceCount / numberOfKernelLaunches2;
  while(sequenceCount > rowsPrKernel2*numberOfKernelLaunches2){
    rowsPrKernel2++;
  }
  printf("rows pr kernel: %d\n",rowsPrKernel2);
  resultMemSize2_dna = rowsPrKernel2 * sequenceCount * sizeof(unsigned int);

  // allocate device memory
  cutilSafeCall(cudaMalloc((void**) &bitStrings2_dna_gpu, bitStringsMemSize2_dna));
  cutilSafeCall(cudaMalloc((void**) &gapFilters2_dna_gpu, bitStringsMemSize2_dna));
  cutilSafeCall(cudaMalloc((void**) &results_dna_gpu, resultMemSize2_dna));
}
