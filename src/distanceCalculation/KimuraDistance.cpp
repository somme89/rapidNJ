#include "KimuraDistance.hpp"
#include "simpleDistanceCalculator.hpp"
#include <math.h>
#include <fstream>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <errno.h>
#include "bitDistanceGap.hpp"
#include "bitDistanceProtein.hpp"
#include "dnaBitString.hpp"

#ifdef __linux__
#include <sys/time.h>
#endif

#ifdef ENABLEGPU
extern "C" void computeDistancesDNA_gpu(unsigned int* ts, unsigned int* tv, unsigned int* gaps);
extern "C" void storeDataDNA_gpu(unsigned int* bitStrings, unsigned int* gapFilters);
extern "C" void initialiseDNA_gpu(unsigned int sequenceCount, unsigned int bitStringCount, unsigned int _bsStride);
extern "C" void computeDistancesProtein_gpu(unsigned int* bitStrings,unsigned int* dist, unsigned int* gaps);
extern "C" void getResultsProtein_gpu(unsigned int* dist, unsigned int* gaps);
extern "C" void initialiseProtein_gpu(unsigned int sequenceCount, unsigned int _bsStride, unsigned int bitStringCount);

extern "C" void computeDistancesDNA2_gpu(float* results);
extern "C" void storeDataDNA2_gpu(unsigned int* bitStrings, unsigned int* gapFilters);
extern "C" void initialiseDNA2_gpu(unsigned int sequenceCount, unsigned int bitStringCount, unsigned int _bsStride);
#endif

static int dayhoff_pams[]={
  195,    196,    197,    198,    199,    200,    200,    201,    202,  203,    
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,    
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,    
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,    
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,    
  246,    248,    249,    250,    252,    253,    254,    255,    257,  258,    
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,    
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,    
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,    
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,    
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,    
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,    
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,    
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,    
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,    
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,    
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,    
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
  988};

KimuraDistance::KimuraDistance(bool verbose, bool fastdist, dataloader* loader, diskMatrix* dm) {
  KimuraDistance::verbose = verbose;
  KimuraDistance::loader = loader;
  KimuraDistance::fastdist = fastdist;
  KimuraDistance::popcnt = popcnt;
  KimuraDistance::seqCount = loader->getSequenceCount();
  KimuraDistance::sequenceNames = *loader->getSequenceNames();  
  KimuraDistance::gpuInitialised = false;
  KimuraDistance::dm = dm;
  distMatrix = NULL;
}

KimuraDistance::KimuraDistance(bool verbose, bool fastdist, dataloader* loader, distType** _distMatrix, diskMatrix* dm) {
  KimuraDistance::verbose = verbose;
  KimuraDistance::loader = loader;
  KimuraDistance::fastdist = fastdist;
  KimuraDistance::popcnt = popcnt;
  KimuraDistance::seqCount = loader->getSequenceCount();
  KimuraDistance::sequenceNames = *loader->getSequenceNames();  
  KimuraDistance::gpuInitialised = false;
  KimuraDistance::dm = dm;
  distMatrix = _distMatrix;
}

void KimuraDistance::computeDistances(int numThreads) {  
  maxDistance = 0.0;
  if(distMatrix != NULL && dm != NULL) {
    cerr << "ERROR: computation of bootstrap values is not supported for trees which use I/O efficient computation." << endl;
    exit(1);
  }
  if(dm == NULL && distMatrix == NULL) {
    distMatrix = new distType*[seqCount];
    for (int i = 0; i < seqCount; i++) {
      distMatrix[i] = new distType[seqCount];
    }
  } else if(dm != NULL) {
    distMatrix = new distType*[numThreads* THREAD_ROW_BUFFER_COUNT];
    for (int i = 0; i < numThreads* THREAD_ROW_BUFFER_COUNT; i++) {
      distMatrix[i] = new distType[seqCount];
    }
  }
  computeDistanceMatrixMT(numThreads);
  if(dm != NULL) {
    for (int i = 0; i < numThreads* THREAD_ROW_BUFFER_COUNT; i++) {
      delete[] distMatrix[i];
    }
    delete[] distMatrix;
  }
}

void KimuraDistance::computeDistancesGPU() { 
#ifndef ENABLEGPU
  cerr << "Computation of distance estimators using GPU is not enabled!" << endl;
  exit(1);
#endif
  if(loader->type == PROTEIN){
    computeDistancesProteinGPU();
  } else {
    computeDistancesDNAGPU2();
  }
}

void KimuraDistance::computeDistancesProteinGPU() { 
#ifdef ENABLEGPU
  unsigned int* bitStrings = loader->bitStringsGPU;
  long paddedLength = (loader->getSequenceLength() + 127) & ~127;

  //allocate memory for the results
  size_t resultMemSize = seqCount * seqCount * sizeof(unsigned int);
  unsigned int* resultDist = (unsigned int*) malloc(resultMemSize);
  unsigned int* resultGap = (unsigned int*) malloc(resultMemSize);

  if(!gpuInitialised){
    initialiseProtein_gpu(seqCount, loader->bsStride, paddedLength);
    gpuInitialised = true;
  }

  computeDistancesProtein_gpu(bitStrings,resultDist,resultGap);
  if(verbose) cerr << "Processing results matrix..." << endl;
  for(int i = 0; i < seqCount-1; i++){
    for (int j = i+1; j < seqCount; j++) {
      distType total = (long) resultDist[i*seqCount + j];
      distType lengthWithoutGaps = paddedLength - resultGap[i*seqCount + j];	      
      distType distance;
      if(lengthWithoutGaps == 0) {
        distance = 0.0f;
      } else {
        distance = total / lengthWithoutGaps;
      }
      if ( distance < 0.75) {
        if (distance > 0.0){
          distance = - log(1.0 - distance - (distance * distance * 0.20));
        }
      } else {
        if (distance > 0.930) {
          distance = 10.0;
        } else {
          int table_index = (int) ((distance*1000.0) - 750.0);
          distance = (distType) (dayhoff_pams[ table_index ] / 100.0);
        }
      }
      if(distance > maxDistance){
        maxDistance = distance;
      }
      distMatrix[i][j] = distance;
    }
  }
  cerr << "Post processing matrix." << endl;
  postProcessDistanceMatrix();
  cerr << "Done!" << endl;
#endif
}

void KimuraDistance::computeDistancesDNAGPU2() {
#ifdef ENABLEGPU
  unsigned int* bitStrings = loader->bitStringsGPU;
  unsigned int* gapFilters = loader->gapFiltersGPU;

  if(!gpuInitialised){
    initialiseDNA2_gpu(seqCount, loader->getBitStringsCount(), loader->bsStride);
    gpuInitialised = true;
  }

  storeDataDNA2_gpu(bitStrings, gapFilters);

  size_t resultMemSize = seqCount * seqCount * sizeof(unsigned int);
  float* results = (float*) malloc(resultMemSize);

  computeDistancesDNA2_gpu(results);

  for(int i = 0; i < seqCount-1; i++){
    for (int j = i+1; j < seqCount; j++) {
      distType distance = results[i*seqCount + j];      
      if(distance > maxDistance){
        maxDistance = distance;
      }
      distMatrix[i][j] = distance;
    }
  }
  postProcessDistanceMatrix();

#endif
}

DistanceEstimate* KimuraDistance::getDistanceEstimateInstance(dataloader* loader) {
  DistanceEstimate* de = NULL;
  if(loader->fastdist) {
    if(loader->type == DNA) {
      de = new bitDistanceGap(loader);
    } else if(loader->type == PROTEIN) {
      de = new bitDistanceProtein(loader);
    } else {
      cerr << "ERROR: Unknown sequence type \"" << loader->type << "\"" << endl;
      exit(1);
    }
  } else {
    de = new simpleDistanceCalculator(loader);
  }
  return de;
}

void KimuraDistance::computeDistanceMatrixMT(int numThreads) { 
  /*  
  pthread_attr_t thread_attr;
  struct sched_param thread_param;
  int thread_policy;
  int ret;
  // set thread priorities CAN ONLY BE EXECUTED AS SUPERUSER 
  pthread_attr_init (&thread_attr);
  ret = pthread_attr_setschedpolicy (&thread_attr, SCHED_FIFO);
  if(ret != 0){
  cout << "ERROR" << endl;
  exit(1);
  }
  ret=  pthread_attr_setinheritsched (&thread_attr, PTHREAD_EXPLICIT_SCHED);
  if(ret != 0){
  cout << "ERROR" << endl;
  exit(1);
  }
  int newPrio = sched_get_priority_max(SCHED_FIFO);
  thread_param.sched_priority = newPrio;
  ret = pthread_attr_setschedparam (&thread_attr, &thread_param);
  if(ret != 0){
  cout << "ERROR" << endl;
  exit(1);
  }
  ret = pthread_setschedparam(pthread_self(),SCHED_FIFO, &thread_param);
  if(ret != 0){
  cout << "ERROR could not set new priority: " << ret << endl;
  cout << strerror(errno) << endl;
  exit(1);
  }
  */

  // Start threads
  pthread_t* threads = new pthread_t[numThreads];
  threadStateKimura** threadStates = new threadStateKimura*[numThreads];
  //start threads
  for (int i = 0; i < numThreads; i++) {
    threadStates[i] = new threadStateKimura();
    threadStates[i]->seqIdx = i;
    threadStates[i]->seqCount = seqCount;
    threadStates[i]->loader = loader;
    threadStates[i]->distMatrix = distMatrix;
    threadStates[i]->numThreads = numThreads;
    threadStates[i]->distanceCalculator = getDistanceEstimateInstance(loader);
    threadStates[i]->availableBuffers = THREAD_ROW_BUFFER_COUNT;
    threadStates[i]->dm = dm;
    if(pthread_mutex_init(&threadStates[i]->mutex, NULL)) {
    cerr << "Could not create mutex" << endl;
    exit(1);
  }
   
    pthread_create(&threads[i], NULL, KimuraDistance::distThread, (void*)threadStates[i]);
  }

  if(dm != NULL) {    
    for(int i = 0; i < seqCount; i++) {
      int threadIdx = i % numThreads;
      while(threadStates[threadIdx]->availableBuffers == THREAD_ROW_BUFFER_COUNT) {
        //busy wait
      }
      //cerr << threadIdx << " " << i << ": " << threadStates[threadIdx]->availableBuffers << endl;
      pthread_mutex_lock(&threadStates[threadIdx]->mutex);
      int bufIdx = i % (numThreads * THREAD_ROW_BUFFER_COUNT);
      //write data
      dm->writeArray(distMatrix[bufIdx], i, seqCount);      
      threadStates[threadIdx]->availableBuffers++;
      pthread_mutex_unlock(&threadStates[threadIdx]->mutex);
    }   
  }
  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
    if(maxDistance < threadStates[i]->maxDistance){
      maxDistance = threadStates[i]->maxDistance;
    }    
  }
  postProcessDistanceMatrix();
  
  for (int i = 0; i < numThreads; i++) {
    pthread_mutex_destroy(&threadStates[i]->mutex);
    delete threadStates[i]->distanceCalculator;
    delete threadStates[i];    
  }
  delete[] threadStates;
  delete[] threads;  
}

// thread entry point
void* KimuraDistance::distThread(void* ptr) {
  distType maxDistance = 0.0;
  threadStateKimura* state = (threadStateKimura*) ptr;
  unsigned long long data[3];
  distType** distMatrix = state->distMatrix;
  int numThreads = state->numThreads;
  unsigned int offset = 0;
  
  for (unsigned int i = state->seqIdx; i < state->seqCount; i+=numThreads) {
    int bufIdx = i;
    if(state->dm != NULL) {
      while(state->availableBuffers == 0) {
        //Busy wait
      }
      bufIdx = i % (numThreads * THREAD_ROW_BUFFER_COUNT);
    } else {
      offset = i+1;
    }

    for (unsigned int j = offset; j < state->seqCount; j++) {
      state->distanceCalculator->computeDistance(i,j,data);
      distType distance = -1.0;
      if(state->loader->type == PROTEIN) {
        long total = data[0] + data[1];	
        if(data[2] != 0) {
          distance = (total / (distType) data[2]);
        }
        if (distance < 0.75) {
          if (distance > 0.0){			
            distance = -log( 1.0 - distance - (distance * distance * 0.20));
          }
        } else {
          if (distance > 0.930) {
            distance = 10.0;
          } else {
            int table_index = (int) ((distance*1000.0) - 750.0);
            distance = (distType) (dayhoff_pams[ table_index ] / 100.0);
          }
        }              
      } else {
        if(data[2] != 0) {
          distType ts = double(data[0]) / double(data[2]);	
          distType tv = double(data[1]) / double(data[2]);
          distType temp1 = 1.0-2.0*ts-tv;
          distType temp2 = 1.0-2.0*tv;
          if(!(temp1 <= 0 || temp2 <= 0)){
            distance = -0.5*log(1.0-2.0*ts-tv) - 0.25*log(1.0-2.0*tv);
          }		  
        }
      }
      if(distance > maxDistance){
        maxDistance = distance;
      }
      distMatrix[bufIdx][j] = distance;
    }
    if(state->dm != NULL) {
      pthread_mutex_lock(&state->mutex);
      state->availableBuffers--;
      pthread_mutex_unlock(&state->mutex);
    }

  }
  state->maxDistance = maxDistance;  
  return NULL;
}

void KimuraDistance::postProcessDistanceMatrix(){
  distType* buffer = NULL;
  if(dm != NULL) {
    buffer = new distType[seqCount];
  }
  for (int i = 0; i < seqCount; i++) {
    if(dm == NULL) {
      buffer = distMatrix[i];  
    } else {
      dm->readArray(buffer, i, seqCount);
    }
    for (int j = 0; j < seqCount; j++) {
      if (j < i && dm == NULL) {
        distMatrix[i][j] = distMatrix[j][i];
      }
      if (i == j) {
        buffer[j] = 0.0;
      }
      if(buffer[j] == -1){
        buffer[j] = maxDistance * 2;
      }
    }
    if(dm != NULL) {
      dm->writeArray(buffer, i, seqCount);
    }
  }
  if(dm != NULL) {
    delete[] buffer;
  }
}

distType** KimuraDistance::getDistanceMatrix(){  
  return distMatrix;
}

KimuraDistance::~KimuraDistance(void){
}
