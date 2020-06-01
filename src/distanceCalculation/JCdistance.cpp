#include "JCdistance.hpp"
#include <math.h>
#include <fstream>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include "simpleDistanceCalculator.hpp"
#include "bitDistanceGap.hpp"
#include "bitDistanceProtein.hpp"

using namespace std;

JCdistance::JCdistance(bool verbose, bool fastdist, dataloader* loader, diskMatrix* dm) {
  JCdistance::verbose = verbose;
  JCdistance::loader = loader;
  JCdistance::fastdist = fastdist;
  JCdistance::seqCount = loader->getSequenceCount();
  JCdistance::seqLength = loader->getSequenceLength();
  JCdistance::sequenceNames = *loader->getSequenceNames();
  JCdistance::dm = dm;
  maxDistance = 0.0; 
}

void JCdistance::computeDistanceMatrix(int numThreads) {  
  maxDistance = 0.0;
  if(dm == NULL) {
    jcDistMatrix = new distType*[seqCount];
    for (int i = 0; i < seqCount; i++) {
      jcDistMatrix[i] = new distType[seqCount];
    }
  } else {    
    jcDistMatrix = new distType*[numThreads* THREAD_ROW_BUFFER_COUNT];
    for (int i = 0; i < numThreads* THREAD_ROW_BUFFER_COUNT; i++) {
      jcDistMatrix[i] = new distType[seqCount];
    }    
  }
  computeDistanceMatrixMT(numThreads);
  if(dm != NULL) {
    for (int i = 0; i < numThreads* THREAD_ROW_BUFFER_COUNT; i++) {
      delete[] jcDistMatrix[i];
    }
    delete[] jcDistMatrix;
  }
}


/*  
pthread_attr_t thread_attr;
struct sched_param thread_param;
int thread_policy;
int ret;

// set thread priorities
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
int newPrio = 99;
thread_param.sched_priority = newPrio;
ret = pthread_attr_setschedparam (&thread_attr, &thread_param);
if(ret != 0){
cout << "ERROR" << endl;
exit(1);
}  
pthread_setschedparam(pthread_self(),SCHED_FIFO, &thread_param);
*/

DistanceEstimate* JCdistance::getDistanceEstimateInstance(dataloader* loader) {
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

void JCdistance::computeDistanceMatrixMT(int numThreads) {   
  // Start threads
  pthread_t* threads = new pthread_t[numThreads];
  threadStateJC** threadStates = new threadStateJC*[numThreads];
  //start threads
  for (int i = 0; i < numThreads; i++) {
    threadStates[i] = new threadStateJC();
    threadStates[i]->seqIdx = i;
    threadStates[i]->seqCount = seqCount;
    threadStates[i]->loader = loader;
    threadStates[i]->jcDistMatrix = jcDistMatrix;
    threadStates[i]->seqLength = seqLength;
    threadStates[i]->numThreads = numThreads;
    threadStates[i]->dm = dm;
    threadStates[i]->availableBuffers = THREAD_ROW_BUFFER_COUNT;
    threadStates[i]->estimator = getDistanceEstimateInstance(loader); 
    if(pthread_mutex_init(&threadStates[i]->mutex, NULL)) {
      cerr << "Could not create mutex" << endl;
      exit(1);
  }    
    //ret = pthread_create(&threads[i], &thread_attr, JCdistance::distJCThread, (void*)threadStates[i]);
    pthread_create(&threads[i], NULL, JCdistance::distJCThread, (void*)threadStates[i]);
  }
  if(dm != NULL) {
    for(int i = 0; i < seqCount; i++) {
      int threadIdx = i % numThreads;
      while(threadStates[threadIdx]->availableBuffers == THREAD_ROW_BUFFER_COUNT) {
        //busy wait
      }      
      pthread_mutex_lock(&threadStates[threadIdx]->mutex);
      int bufIdx = i % (numThreads * THREAD_ROW_BUFFER_COUNT);
      //write data
      dm->writeArray(jcDistMatrix[bufIdx], i, seqCount);
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

  for (int i = 0; i < numThreads; i++) {
    pthread_mutex_destroy(&threadStates[i]->mutex);
    delete threadStates[i]->estimator;
    delete threadStates[i];    
  }
  delete[] threadStates;
  delete[] threads;
  
  postProcessDistanceMatrix();
}

//thread entry point
void* JCdistance::distJCThread(void* ptr) {
  distType maxDistance = 0.0;
  threadStateJC* state = (threadStateJC*) ptr;
  unsigned long long distances[3];
  distType** jcDistMatrix = state->jcDistMatrix;
  int numThreads = state->numThreads;
  DistanceEstimate* estimator = state->estimator;
  double alpha = 4.0;
  if(state->loader->type == PROTEIN) {
	alpha = 20.0;
  }
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
      estimator->computeDistance(i,j,distances);
      long total = distances[0] + distances[1];
      distType distance;
      if(distances[2] == 0){
        distance = 0;
      } else {
        distance = total / (distType) distances[2];
      }
	  if(distance / ((alpha-1)/alpha) >= 1.0) {
		distance = -1;		
	  } else {		  
		 distance = -((alpha-1)/alpha) * log( 1.0 - distance / ((alpha-1)/alpha));
	  }
      if(maxDistance < distance){
        maxDistance = distance;
      }
      jcDistMatrix[bufIdx][j] = distance;	  
    }	
    if(state->dm != NULL) {
      pthread_mutex_lock(&state->mutex);
      state->availableBuffers--;
      pthread_mutex_unlock(&state->mutex);
    }
  }
  state->maxDistance = maxDistance;
  delete estimator;
  state->estimator = NULL;
  return NULL;  
}

void JCdistance::postProcessDistanceMatrix(){
  distType* buffer = NULL;
  if(dm != NULL) {
    buffer = new distType[seqCount];
  }
  for (int i = 0; i < seqCount; i++) {
    if(dm == NULL) {
      buffer = jcDistMatrix[i];  
    } else {
      dm->readArray(buffer, i, seqCount);
    }
    for (int j = 0; j < seqCount; j++) {
      if (j < i && dm == NULL) {
        jcDistMatrix[i][j] = jcDistMatrix[j][i];
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

distType** JCdistance::getDistanceMatrix(){  
  return jcDistMatrix;
}
