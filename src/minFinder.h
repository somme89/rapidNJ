#ifndef MINFINDER_H_
#define MINFINDER_H_
#include "threadedNJ.h"

class minFinder
{
	
public:
	minFinder(threadState* state);
	void run();
	void findMin();
	~minFinder();
	
private:
	threadState* state;
	
};

#endif /*MINFINDER_H_*/
