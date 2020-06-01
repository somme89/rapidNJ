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

__device__ void parallelReductionDNA(){
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
__global__ void computeSliceDNA(unsigned int* gaps_gpu, unsigned int* ts_gpu, unsigned int* tv_gpu, unsigned int* bitStrings_gpu, unsigned int* gapFilters_gpu, unsigned int bsStride, unsigned int dataSize, unsigned int rowOffset){

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
  parallelReductionDNA();

  //store data in global memory
  if(threadIdx.x == 0){
    //printf("%d \n", blockIdx.x * sequenceCount_GPU + i);
    tv_gpu[blockIdx.x * gridDim.y + i] = shared_data_dna[0];
    ts_gpu[blockIdx.x * gridDim.y + i] = shared_data_dna[threads_pr_block];
    gaps_gpu[blockIdx.x * gridDim.y + i] = shared_data_dna[threads_pr_block*2];
  }
}

unsigned int sequenceCount_dna;
unsigned int dataSize_dna;
size_t resultMemSize_dna;
size_t bitStringsMemSize_dna;
//results
unsigned int* ts_gpu;
unsigned int* tv_gpu;
unsigned int* gaps_dna_gpu;
unsigned int bsStride_dna;
//data
unsigned int* gapFilters_dna_gpu;
unsigned int* bitStrings_dna_gpu;
unsigned int rowsPrKernel;
unsigned int numberOfKernelLaunches;

#ifdef TIMING
extern float totalGpuComputation;
extern timeval start,end;
extern float totalTransfer;
#endif

extern "C" void computeDistancesDNA_gpu(unsigned int* ts, unsigned int* tv, unsigned int* gaps) {

  // execute the kernel
  //  unsigned int startIdx = i * threads_pr_block;
  for(int i = 0; i < numberOfKernelLaunches; i++){
    printf("kernel: %d -----------------------\n",i);
    // setup execution parameters
    unsigned int gridSize_dna = 0;
    if(i != numberOfKernelLaunches-1){
      gridSize_dna = rowsPrKernel;
    } else {
      gridSize_dna = sequenceCount_dna - (numberOfKernelLaunches-1)*rowsPrKernel;
    }
    dim3 grid(gridSize_dna,sequenceCount_dna);
    dim3 block(threads_pr_block, 1, 1);
    unsigned int transferSize = gridSize_dna * sequenceCount_dna * sizeof(unsigned int);
    computeSliceDNA<<< grid, block >>>(gaps_dna_gpu, ts_gpu, tv_gpu, bitStrings_dna_gpu, gapFilters_dna_gpu, bsStride_dna, dataSize_dna, i*rowsPrKernel);
    
#ifdef TIMING
    cudaThreadSynchronize();
    gettimeofday(&end,NULL);
    totalGpuComputation += (end.tv_sec - start.tv_sec)*1000.0 + (end.tv_usec - start.tv_usec)/1000.0;
    gettimeofday(&start,NULL);
#endif

  cutilSafeCall(cudaMemcpy(&ts[i * rowsPrKernel * sequenceCount_dna], ts_gpu, transferSize, cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaMemcpy(&tv[i * rowsPrKernel * sequenceCount_dna], tv_gpu, transferSize, cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaMemcpy(&gaps[i * rowsPrKernel * sequenceCount_dna], gaps_dna_gpu, transferSize, cudaMemcpyDeviceToHost));
  
#ifdef TIMING
  cudaThreadSynchronize();
  gettimeofday(&end,NULL);
  totalTransfer += (end.tv_sec - start.tv_sec)*1000.0 + (end.tv_usec - start.tv_usec)/1000.0;
  gettimeofday(&start,NULL);
#endif
  }
}

extern "C" void getResultsDNA_gpu(unsigned int* ts, unsigned int* tv, unsigned int* gaps) {
  // Wait for kernels to finish.
  cutilCheckMsg("Kernel execution failed");  

  // Copy results from device to host
  /*cutilSafeCall(cudaMemcpy(ts, ts_gpu, resultMemSize_dna, cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaMemcpy(tv, tv_gpu, resultMemSize_dna, cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaMemcpy(gaps, gaps_dna_gpu, resultMemSize_dna, cudaMemcpyDeviceToHost));*/
}

extern "C" void storeDataDNA_gpu(unsigned int* bitStrings, unsigned int* gapFilters) {
  // Copy results from device to host
  cutilSafeCall(cudaMemcpyAsync(bitStrings_dna_gpu, bitStrings, bitStringsMemSize_dna, cudaMemcpyHostToDevice,0));
  cutilSafeCall(cudaMemcpyAsync(gapFilters_dna_gpu, gapFilters, bitStringsMemSize_dna, cudaMemcpyHostToDevice,0));
}

extern "C" void initialiseDNA_gpu(unsigned int sequenceCount, unsigned int bitStringCount, unsigned int _bsStride) {
  printf("initialising GPU... \n");
  sequenceCount_dna = sequenceCount;
  bsStride_dna = _bsStride;  
  dataSize_dna = min(bsStride_dna,bitStringCount * 4);

  //printf("TEST: %d %d %d\n",bsStride_dna, bitStringCount, dataSize_dna);
  //exit(0);
  
  
  resultMemSize_dna = sequenceCount * sequenceCount * sizeof(unsigned int);
  bitStringsMemSize_dna = bsStride_dna * sizeof(unsigned int) * sequenceCount;
  
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);  
  float freeGpuMem = deviceProp.totalGlobalMem - 250*1024*1024; //substract 100MB
  
    
  //printf("TODO REMOVE THIS \n");
  //float freeGpuMem = 50.0 * 1024 * 1024;

  freeGpuMem -= bitStringsMemSize_dna;
  float temp = min(1.0,freeGpuMem / (resultMemSize_dna*3.0));
  rowsPrKernel = sequenceCount * temp;
  numberOfKernelLaunches = sequenceCount/rowsPrKernel;
  if(sequenceCount % rowsPrKernel != 0){
    numberOfKernelLaunches++;
  }

  rowsPrKernel = sequenceCount / numberOfKernelLaunches;
  unsigned int exRows = sequenceCount - rowsPrKernel*numberOfKernelLaunches;
  rowsPrKernel += exRows;

  printf("rows pr kernel: %d\n",rowsPrKernel);
  resultMemSize_dna = rowsPrKernel * sequenceCount * sizeof(unsigned int);

  //exit(0);
  // allocate device memory
  cutilSafeCall(cudaMalloc((void**) &bitStrings_dna_gpu, bitStringsMemSize_dna));
  cutilSafeCall(cudaMalloc((void**) &gapFilters_dna_gpu, bitStringsMemSize_dna));
  cutilSafeCall(cudaMalloc((void**) &ts_gpu, resultMemSize_dna));
  cutilSafeCall(cudaMalloc((void**) &tv_gpu, resultMemSize_dna));
  cutilSafeCall(cudaMalloc((void**) &gaps_dna_gpu, resultMemSize_dna));
}
