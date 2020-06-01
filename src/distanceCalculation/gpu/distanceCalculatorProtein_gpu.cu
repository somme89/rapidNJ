#include <cutil_inline.h>
#include "constants.hpp"
#include <math.h>
#include <sys/time.h>

__shared__ unsigned int shared_data_protein[threads_pr_block];

__device__ void parallelReductionProtein(unsigned int val){
  //wait for all threads to finish and use shared memory to store results
  shared_data_protein[threadIdx.x] = val;
  __syncthreads();

   //perform parallel reduction
  for(unsigned int s=blockDim.x/2; s>32; s>>=1){    
    if (threadIdx.x < s){
      shared_data_protein[threadIdx.x ] += shared_data_protein[threadIdx.x + s];
    } 
    __syncthreads();
  }
  if (threadIdx.x < 32) {
    shared_data_protein[threadIdx.x ] += shared_data_protein[threadIdx.x + 32 ];
  }   
  __syncthreads();
  

  if (threadIdx.x < 16) {
    shared_data_protein[threadIdx.x ] += shared_data_protein[threadIdx.x + 16 ];
  }
  __syncthreads();
  if (threadIdx.x < 8) {
    shared_data_protein[threadIdx.x ] += shared_data_protein[threadIdx.x +  8 ];
  }
  __syncthreads();

  if (threadIdx.x < 4) {
    shared_data_protein[threadIdx.x ] += shared_data_protein[threadIdx.x +  4 ];
  }  
  __syncthreads();
  
  if (threadIdx.x < 2) {
    shared_data_protein[threadIdx.x ] += shared_data_protein[threadIdx.x + 2 ];
  }  
  __syncthreads();

  if(threadIdx.x == 0){
    shared_data_protein[0 ] += shared_data_protein[1 ];
  }
  __syncthreads();
}

__global__ void computeDistanceProtein(unsigned int* gaps_gpu, unsigned int* dist_gpu, unsigned int* bitStrings_gpu, unsigned int bsStride, unsigned int rowOffset, unsigned int dataSize_protein){
    
  int seq2 = blockIdx.y;
  if(seq2 < blockIdx.x+1+rowOffset){
    return;
  }

  uchar4* charData = (uchar4*) bitStrings_gpu;
  unsigned int iterations = dataSize_protein / threads_pr_block;
  if(dataSize_protein % threads_pr_block != 0){
    iterations++;
  }
  //printf(" %d \n",iterations);
  //for(int seq2 = blockIdx.x+1+rowOffset; seq2 < sequenceCount_gpu; seq2++) {
  unsigned int sumR = 0;
  unsigned int sumG = 0;
  for(int slice = 0; slice < iterations; slice++) {      
    uchar4 r, g;
    // fetch next block of sequence data
    unsigned int idx0 = (rowOffset+blockIdx.x)*bsStride + slice*threads_pr_block + threadIdx.x;
    uchar4 base = charData[idx0];
    unsigned int idx1 = seq2*bsStride + slice*threads_pr_block + threadIdx.x;
    uchar4 target = charData[idx1];
    
    r.x = base.x != target.x;
    g.x = base.x < 64 || target.x < 64;
    r.x = r.x && !g.x;
    
    r.y = base.y != target.y;
    g.y = base.y < 64 || target.y < 64;
    r.y = r.y && !g.y;
    
    r.z = base.z != target.z;
    g.z = base.z < 64 || target.z < 64;
    r.z = r.z && !g.z;
    
    r.w = base.w != target.w;
    g.w = base.w < 64 || target.w < 64;
    r.w = r.w && !g.w;
    
    if(slice * threads_pr_block + threadIdx.x < dataSize_protein) {
      sumR += r.x;
      sumR += r.y;
      sumR += r.z;
      sumR += r.w;
      
      sumG += g.x;
      sumG += g.y;
      sumG += g.z;
      sumG += g.w;
      //printf("%d: %d %d: %d %d %d %d\n",bsStride,idx0,idx1, g.x, g.y, g.z, g.w);
    }
  }
  //sum over all threads
  parallelReductionProtein(sumR);
  sumR = shared_data_protein[0];
  parallelReductionProtein(sumG);
  sumG = shared_data_protein[0];
  
  //store result in global memory
  if(threadIdx.x == 0) {
    dist_gpu[blockIdx.x * gridDim.y + seq2] = sumR;
    gaps_gpu[blockIdx.x * gridDim.y + seq2] = sumG;
  }
}

//-------------------------------------------------------------------------------------------------------------

unsigned int gridSize_protein;
size_t resultMemSize_protein;
size_t bitStringsMemSize_protein;
unsigned int sequenceCount;
//results
unsigned int* dist_protein_gpu;
unsigned int* gaps_protein_gpu;
unsigned int bsStride_protein;
//data
unsigned int* bitStrings_protein_gpu;
unsigned int numberOfKernelLaunches_protein;
unsigned int rowsPrKernel_protein;
unsigned int dataSize_protein;

#ifdef TIMING
extern float totalGpuComputation;
extern timeval start,end;
extern float totalTransfer;
#endif

extern "C" void computeDistancesProtein_gpu(unsigned int* bitStrings, unsigned int* dist, unsigned int* gaps) {
  //  printf("%p %p %d \n",bitStrings_protein_gpu,bitStrings,bitStringsMemSize_protein);
  cutilSafeCall(cudaMemcpyAsync(bitStrings_protein_gpu, bitStrings, bitStringsMemSize_protein, cudaMemcpyHostToDevice,0));
  for(int i = 0; i < numberOfKernelLaunches_protein; i++){
    printf("kernel: %d -----------------------\n",i);
    // setup execution parameters
    unsigned int gridSize = 0;
    if(i != numberOfKernelLaunches_protein-1){
      gridSize = rowsPrKernel_protein;
    } else {
      gridSize = sequenceCount - (numberOfKernelLaunches_protein-1)*rowsPrKernel_protein;
    }
    printf("gridSize: %d %d\n",gridSize,sequenceCount);
    dim3 grid(gridSize,sequenceCount);
    dim3 block(threads_pr_block, 1, 1);
    unsigned int transferSize = gridSize * sequenceCount * sizeof(unsigned int);
    computeDistanceProtein<<< grid, block >>>(gaps_protein_gpu, dist_protein_gpu, bitStrings_protein_gpu, bsStride_protein, i*rowsPrKernel_protein, dataSize_protein);

#ifdef TIMING
    cudaThreadSynchronize();
    gettimeofday(&end,NULL);
    totalGpuComputation += (end.tv_sec - start.tv_sec)*1000.0 + (end.tv_usec - start.tv_usec)/1000.0;
    gettimeofday(&start,NULL);
#endif

    printf("Copying results...\n");
    cutilSafeCall(cudaMemcpy(&dist[i * rowsPrKernel_protein * sequenceCount], dist_protein_gpu, transferSize, cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(&gaps[i * rowsPrKernel_protein * sequenceCount], gaps_protein_gpu, transferSize, cudaMemcpyDeviceToHost));
    printf("Finished\n");

#ifdef TIMING
    cudaThreadSynchronize();
    gettimeofday(&end,NULL);
    totalTransfer += (end.tv_sec - start.tv_sec)*1000.0 + (end.tv_usec - start.tv_usec)/1000.0;
    gettimeofday(&start,NULL);
#endif
  }
}

extern "C" void getResultsProtein_gpu(unsigned int* dist, unsigned int* gaps) {
  /*  // Wait for kernels to finish.
  cutilCheckMsg("Kernel execution failed");  
  
  // Copy results from device to host
  cutilSafeCall(cudaMemcpy(dist, dist_protein_gpu, resultMemSize_protein, cudaMemcpyDeviceToHost));
  cutilSafeCall(cudaMemcpy(gaps, gaps_protein_gpu, resultMemSize_protein, cudaMemcpyDeviceToHost));
  */
  printf("getResult is not implemented for proteins \n");
  exit(1);
}

extern "C" void initialiseProtein_gpu(unsigned int _sequenceCount, unsigned int _bsStride, unsigned int paddedLength) {
  //  printf("initialising GPU... \n");
  bsStride_protein = _bsStride;
  sequenceCount = _sequenceCount;
  dataSize_protein = paddedLength / 4;  
  //printf("%d %d \n",bsStride_protein,paddedLength);
  //exit(0);

  resultMemSize_protein = sequenceCount * sequenceCount * sizeof(unsigned int);
  bitStringsMemSize_protein = bsStride_protein * sizeof(unsigned int) * sequenceCount;
  
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);  
  float freeGpuMem = deviceProp.totalGlobalMem - 250*1024*1024; //substract 100MB
    
  freeGpuMem -= bitStringsMemSize_protein;
  float temp = min(1.0,freeGpuMem / (resultMemSize_protein*2.0));
  rowsPrKernel_protein = sequenceCount * temp;

  //rowsPrKernel_protein = 500;
  
  numberOfKernelLaunches_protein = sequenceCount/rowsPrKernel_protein;
  if(sequenceCount % rowsPrKernel_protein != 0){
    numberOfKernelLaunches_protein++;
  }
  
  printf("Number of kernel launches needed: %d\n",numberOfKernelLaunches_protein);
  
  rowsPrKernel_protein = sequenceCount / numberOfKernelLaunches_protein;
  unsigned int exRows = sequenceCount - rowsPrKernel_protein*numberOfKernelLaunches_protein;
  rowsPrKernel_protein += exRows;
  
  printf("rows pr kernel: %d\n",rowsPrKernel_protein);
  resultMemSize_protein = rowsPrKernel_protein * sequenceCount * sizeof(unsigned int);
  
  // allocate device memory
  cutilSafeCall(cudaMalloc((void**) &bitStrings_protein_gpu, bitStringsMemSize_protein)); 
  cutilSafeCall(cudaMalloc((void**) &dist_protein_gpu, resultMemSize_protein));
  cutilSafeCall(cudaMalloc((void**) &gaps_protein_gpu, resultMemSize_protein));
}
