#ifndef __DISTANCE_ESTIMATE_HPP__
#define __DISTANCE_ESTIMATE_HPP__

#include <stdinclude.h>
#include "dataloader.hpp"

using namespace std;

class DistanceEstimate {

public:
  DistanceEstimate(dataloader* loader);
  virtual ~DistanceEstimate();
  virtual void computeDistance(int i, int j, unsigned long long* data) = 0;

protected:
  bool fastDist;
  dataloader* loader;
  InputType type;
};

#endif
