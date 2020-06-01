#ifndef THREADEDNJ_H_
#define THREADEDNJ_H_

#include "node.h"
#include "datareader.h"

class threadedNJ
{

public:
  threadedNJ(datareader* reader);
  node* run();
  ~threadedNJ();
  static void* njThread( void* ptr );    	
  int min1;
  int min2;    
  
  
 private:
  void findMin();
  
  void initialize();
  void mergeMinNodes();
  void updateMatrix();
};

struct threadState {	
  int rowstart;
  int rowEnd;  
  int min1;
  int min2;
  double min;
  bool run;
};

void executefind(threadState* state);

#endif /*THREADEDNJ_H_*/
