#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <map>

using namespace std;

map<int,char> nuc_map;

void simPhylipSequences(int seq_num, int seq_length, bool useGaps){
  srand(time(NULL));
  int total_char = seq_length, remaining_char = seq_length;
  char* rnd_seq = new char[60];
  for(int i= 0; i < 60; i++){
    rnd_seq[i]= nuc_map[rand() % nuc_map.size()];
  }
  cout << ' ' << seq_num << ' ' << seq_length << endl;
  char buffer[16];
  while(remaining_char > 0){
    for(int i = 0; i < seq_num; i++){
      //change one random nucletide in each sequence
      int ranInt = rand() % 60;
      char original_nuc = rnd_seq[ranInt];
      rnd_seq[ranInt] = nuc_map[rand() % nuc_map.size()];
      sprintf(buffer,"seq_%06d ",i);
      cout << buffer << ' ';
      if(remaining_char <=60){
	for(int j = 0; j < remaining_char; j++){
	  cout << rnd_seq[j];
	}
      }	else{
	for(int j = 0; j < 60; j++){
	  cout << rnd_seq[j];
	}
      }
      cout << endl;
      rnd_seq[ranInt] = original_nuc;
    }
    remaining_char -= 60;
    //generate new random sequence
    for(int i= 0; i < 60; i++){
      rnd_seq[i]= nuc_map[rand() % nuc_map.size()];
    }
    cout << endl;
    cerr << (1-remaining_char/ ((float)total_char))*100.0 << '\r';
  }
}

void simStockholmSequences(int seq_num, int seq_length, bool useGaps){
  double changeProb = 30;  
  srand(time(NULL));
  // build random sequence
  char* rnd_seq = new char[seq_length];
  for(int i= 0; i < seq_length; i++){
    rnd_seq[i]= nuc_map[rand() % nuc_map.size()];
  }
  char buffer[16];
  cout << "#STOCKHOLM 1.0" << endl;
  for(int i = 0; i < seq_num; i++){
    sprintf(buffer,"seq_%06d ",i);
    cout << buffer << '\t';    
    for(int j = 0; j < seq_length; j++){
      if(rand() % 100 > changeProb){
	cout << rnd_seq[j];
      } else {
	cout << nuc_map[rand() % nuc_map.size()];
      }
    }
    cout << endl;
  }
  cout << "//" << endl;
}

void simFastaSequences(int seq_num, int seq_length, bool useGaps){
  double changeProb = 30;  
  srand(time(NULL));
  // build random sequence
  char* rnd_seq = new char[seq_length];
  for(int i= 0; i < seq_length; i++){
    rnd_seq[i]= nuc_map[rand() % nuc_map.size()];
  }
  char buffer[16];
  for(int i = 0; i < seq_num; i++){
    sprintf(buffer,">seq_%06d ",i);
    cout << buffer << endl;    
    for(int j = 0; j < seq_length; j++){
      if(j != 0 && j % 60 == 0){
	cout << endl;
      }
      if(rand() % 100 > changeProb){
	cout << rnd_seq[j];
      } else {
	cout << nuc_map[rand() % nuc_map.size()];
      }      
    }
    cout << endl;
  }
}

int main(int argc, char* argv[]){
  if(argc < 5 || argc > 6){
    cout << "USAGE: sequence_count sequence_length type format [-g]" << endl;
    cout << "type: p=protein d=DNA" << endl;
    cout << "format: p=phylip, s=stockholm, f=fasta" << endl;
    exit(0);
  }
  int seq_num = atoi(argv[1]);
  int seq_length = atoi(argv[2]);
  string type = argv[3];
  string format = argv[4];
  bool useGaps = false;  
  if(argc == 6){
    if(string(argv[5]) == "-g"){
      useGaps = true;
    } else {
      cout << "unrecognised argument " << argv[5] << endl;
      exit(0);
    } 
  }
  if(type == "d"){
    nuc_map[0] = 'A';
    nuc_map[1] = 'C';
    nuc_map[2] = 'G';
    nuc_map[3] = 'T';
    if(useGaps){
      nuc_map[4] = '-';
    }
  } else if(type == "p"){
    nuc_map[0] = 'A';
    nuc_map[1] = 'R';
    nuc_map[2] = 'N';
    nuc_map[3] = 'D';
    nuc_map[4] = 'C';
    nuc_map[5] = 'E';
    nuc_map[6] = 'Q';
    nuc_map[7] = 'G';
    nuc_map[8] = 'H';
    nuc_map[9] = 'I';
    nuc_map[10] = 'L';
    nuc_map[11] = 'K';
    nuc_map[12] = 'M';
    nuc_map[13] = 'F';
    nuc_map[14] = 'P';
    nuc_map[15] = 'S';
    nuc_map[16] = 'T';
    nuc_map[17] = 'W';
    nuc_map[18] = 'Y';
    nuc_map[19] = 'V';
    if(useGaps){
      nuc_map[20] = '-';
    }
  } else {
    cout << "Unknown type " << type << endl;
  }
  if(format == "p"){
    simPhylipSequences(seq_num, seq_length, useGaps);
  } else if(format == "s") {
    simStockholmSequences(seq_num, seq_length, useGaps);
  } else if(format == "f"){
    simFastaSequences(seq_num, seq_length, useGaps);
  } else {
    cerr << "Unknown format" << endl;
    exit(1);
  }
}
