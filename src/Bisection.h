#ifndef BISECTION_H
#define BISECTION_H
#include "stdinclude.h"
#include <algorithm>

class bisection {
public:
  vector<int>* section1;
  vector<int>* section2;
  int node1;
  int node2;
  int hashCached;
  int isHashCached;

  bisection(){
    section1 = NULL;
    section2 = NULL;
    isHashCached = 0;
  }

  bisection(int n1, int n2){
    node1 = n1;
    node2 = n2;
    isHashCached = 0;
    section1 = new vector<int>();
    section2 = new vector<int>();
  }

  int hash(){
    if(isHashCached == 1) {
      return hashCached;
    }
    int hash = 0;
    if(section1->size() > section2->size()) {
      for(int i = 0; i < section1->size(); i++) {
        hash += section1->at(i);
      }
    } else if (section1->size() < section2->size()){
      for(int i = 0; i < section2->size(); i++) {
        hash += section2->at(i);
      }
    } else {
      for(int i = 0; i < section1->size(); i++) {
        hash += section1->at(i);
      }
      for(int i = 0; i < section2->size(); i++) {
        hash += section2->at(i);
      }
    }     
    return hash;
  }

  bool compareSections(vector<int>* s1, vector<int>* s2) {
    if(s1->size() != s2->size()) {
      return false;
    }
    sort(s1->begin(), s1->end());
    sort(s2->begin(), s2->end());
    for(int i = 0; i < s1->size(); i++) {
      if(s1->at(i) != s2->at(i)) {
        return false;
      }
    }
    return true;
  }

  bool equals(bisection* other){
    //Only need to find one match as this implies that the other section will also match.
    if(compareSections(section1, other->section1)) {
      return true;
    } else if(compareSections(section1, other->section2)) {
      return true;
    } else {
      return false;
    }
  }

  ~bisection(){
    if(section1 != NULL) {
      delete section1;
    }
    if(section2 != NULL) {
      delete section2;
    }
  }

};

#endif
