#include "../stdinclude.h"
#include "dnaBitString.hpp"

ostream& operator<<(ostream& out, const dnaBitString& s) {
    for (int i = 0;i < 4 ; i++)  {
        cout << "[";
        for (int j = 0;j < 16 ; j++) {
            if(j%2 == 0 && j != 0){
                    cout << "-";
            }
	    out << ((s.data.d[i]  >> (j*2)) & 1);
            out << ((s.data.d[i]  >> (j*2+1)) & 1);
        }
        out << "] ";
    }
    return out;
}

dnaBitString::dnaBitString() {}

/**
* Builds a 128 bit sequence from an amino acid sequence.
* [(part0),(part1),(part2),(part3)]
* each part is 32 bit and represents 16 characters in reverse
* ACGT is represented by the bit string 00000000 00000000 00000000 10011100
*/
dnaBitString::dnaBitString(char* seq) {
    for (int idx = 0; idx < 4; idx++) {
        data.d[idx] = 0;
        for (int i = 0; i < 16; i++) {
            switch (seq[i + (idx*16)]) {
            //A has value 00 so nothing is added
            case 'C': {
                data.d[idx] += (Cbin<<(i*2));
                break;
            }
            case 'G': {
                data.d[idx] += (Gbin<<(i*2));
                break;
            }
            case 'T': {
                data.d[idx] += (Tbin<<(i*2));
                break;
            }
            default:
                break;
            }
        }
    }
}

dnaBitString::dnaBitString(v4ui initData){
  data.v = initData;
}


void dnaBitString::printSequence() {
    for (int i = 0;i < 4 ; i++)  {
        for (int j = 0;j < 16 ; j++) {
            switch ((data.d[i]  >> (2*j) ) & 3) {
            case 0: {
                cout << 'A';
                break;
            }
            case 1: {
                cout << 'G';
                break;
            }
            case 2: {
                cout << 'T';
                break;
            }
            case 3: {
                cout << 'C';
                break;
            }
            default: {
                cout << '?';
                }
            }
        }
        cout << " ";
    }
}




