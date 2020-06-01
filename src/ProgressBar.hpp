#pragma once
#include <stack>
#include <utility>
#include "stdinclude.h"

class ProgressBar
{
public:
  ProgressBar(void);
  void childProgress(double fraction);
  void setProgress(double progress);
  void finish();
  ~ProgressBar(void);
  
private:  
  std::stack<std::pair<double,double> > intervals;
  double curStart;
  double curEnd;
  double currentProgress;
};

