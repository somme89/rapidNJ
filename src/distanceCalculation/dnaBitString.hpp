#ifndef DNABITSTRING_H
#define DNABITSTRING_H

#include <iostream>

using namespace std;

class dnaBitString {

public:
  dnaBitString();
  dnaBitString(char*);
  dnaBitString(v4ui initData);
  void printSequence();
  
  v4ui_vector data;
  
private:
    friend ostream& operator<<(ostream& out, const dnaBitString& s) ;
};

ostream& operator<<(ostream& out, const dnaBitString& s) ;
#endif
