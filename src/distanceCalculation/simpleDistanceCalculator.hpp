#ifndef SIMPLEDISTANCECALCULATOR_H
#define SIMPLEDISTANCECALCULATOR_H

#include "stdinclude.h"
#include "dataloader.hpp"
#include "DistanceEstimate.hpp"

class simpleDistanceCalculator : public DistanceEstimate {

 public:
  simpleDistanceCalculator(dataloader* loader);
  ~simpleDistanceCalculator();
  void computeDistance(int i, int j, unsigned long long* data);
 
 private:
  vector<char*> sequences;
  std::string* sequenceNames;
  int seqLength;
  void calculateTsTv(int i, int j, unsigned long long* data);
  void calculateHammingDistance(int i, int j, unsigned long long* data);
};

#endif
