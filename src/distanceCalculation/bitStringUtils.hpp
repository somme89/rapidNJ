#ifndef BITSTRINGUTILS_HPP
#define BITSTRINGUTILS_HPP
#include "stdinclude.h"

namespace {

  void printBitString(v4ui data){
    v4ui_vector vector;
    vector.v = data;
    for (int i = 0;i < 4 ; i++)  {
      cout << "[";
      for (int j = 0;j < 32 ; j++) {
        if(j%8 == 0 && j != 0){
          cout << " ";
        }
        cout << ((vector.d[i]  >> j) & 1);
      }
      cout << "] ";
    }
    cout << endl;
  }

  void printBitString(unsigned int data){
    cout << "[";
    for (int j = 0;j < 32 ; j++) {
      if(j%8 == 0 && j != 0){
        cout << " ";
      }
      cout << ((data  >> j) & 1);
    }
    cout << "] ";    
    cout << endl;
  }

  void printProteinSequence(v4ui data){
    v4ui_vector vector;
    vector.v = data;
    for (int i = 0;i < 4 ; i++)  {
      cout << "[";
      for (int j = 0;j < 4; j++) {
        if(j != 0){
          cout << " ";
        }
        unsigned int c = ((vector.d[i] >> (j*8)) & 255);
        //cout << (char) c << "       ";
        cout << c << "       ";
      }
      cout << "] ";
    }
    cout << endl;
  }

  void printDNASequence(v4ui data){
    v4ui_vector vector;
    vector.v = data;
    for (int i = 0;i < 4 ; i++)  {
      cout << "[";
      for (int j = 0;j < 16; j++) {
        if(j % 4 == 0 && j != 0) {
          cout << " ";
        }
        unsigned int c = ((vector.d[i] >> (j*2)) & 3);
        switch (c) {
        case 0: {
          cout << "A ";
          break;
                }
        case 1: {
          cout << "G ";
          break;
                }
        case 2: {
          cout << "T ";
          break;
                }
        case 3: {
          cout << "C ";
          break;
                }	  
        default:
          break;
        }
      }
      cout << "] ";
    }
    cout << endl;
  }

  void printDNASequence(unsigned int data){
    cout << "[";
    for (int j = 0;j < 16; j++) {
      if(j % 4 == 0 && j != 0) {
        cout << " ";
      }
      unsigned int c = ((data >> (j*2)) & 3);
      switch (c) {
      case 0: {
        cout << "A ";
        break;
              }
      case 1: {
        cout << "G ";
        break;
              }
      case 2: {
        cout << "T ";
        break;
              }
      case 3: {
        cout << "C ";
        break;
              }	  
      default:
        break;
      }
    }
    cout << "] ";
    cout << endl;
  }

}

#endif
