#include "ProgressBar.hpp"
#include <iomanip>

ProgressBar::ProgressBar(void)
{
  curStart = 0;
  curEnd = 100;  
  currentProgress = 0;
  //intervals.push(make_pair(curStart, curEnd));
}

void ProgressBar::childProgress(double fraction)
{
  intervals.push(make_pair(curStart, curEnd));  
  curEnd = (curEnd - curStart) * fraction + currentProgress;
  curStart = currentProgress;
  if(curStart > 100) {
    curStart = 100;
  }
  if(curEnd > 100) {
    curEnd = 100;
  }
  if(curStart > curEnd) {
    curEnd = curStart;
  }
}

void ProgressBar::setProgress(double progress) {
  if(progress > 1 || progress < 0) {
    cerr << "INTERNAL ERROR: Progress has to be in [0;1]" << endl;
    exit(1);
  }
  double intervalProgress = (curEnd - curStart) * progress;
  currentProgress = curStart + intervalProgress;
  if(currentProgress > curEnd) {
    currentProgress = curEnd;
  }
  cerr << setiosflags(ios::fixed) << setprecision(2) << currentProgress << "% \r";
  cerr.flush();
}

void ProgressBar::finish() {
  setProgress(1);
  if(intervals.size() > 0) {    
    pair<double,double> temp = intervals.top();
    curStart = temp.first;
    curEnd = temp.second;
    intervals.pop();
  }
}

ProgressBar::~ProgressBar(void)
{
}
