/*A simple implementation of the neighbour-joining method*/

#include "stdinclude.h"
#include "testNJ.h"
#include "float.h"

using namespace std;

void printMatrix(distType** matrix, int size);
void printArray(distType* a, int size);

testNJ::testNJ(distMatrixReader* reader, int matrixSize)
{
  testNJ::matrixSize = matrixSize;
  testNJ::reader = reader;
  matrix = reader->getMatrix();
  separationsums = new distType[matrixSize];
  separations = new distType[matrixSize];
  clusterCount = matrixSize;      
  activeRows = new int[matrixSize*2];
  nextId = matrixSize;
  initialize();
}

void testNJ::step(int idx1, int idx2, distType value){
  //  cout << activeRows[38] << endl;
  if( clusterCount == 2){
    // finish by joining the two remaining clusters
    min1 = -1;
    min2 = -1;
    // find the last nodes
    for(int i = 0; i < matrixSize; i++){
      if(activeRows[i] != -1){
	if(min1 == -1){
	  min1 = i;
	} else {
	  min2 = i;	
	  break;
	}            
      }
    }  
    double distance = matrix[min1][min2];
    if(idx1 == min1 && idx2 == min2){
      if(distance != value){
	cout << "ERROR: indexes match but value doesn't: " << value << "!=" << min << endl;  
	exit(1);
      }
    } else if(value != distance){
      cout << "ERROR: both value and indexes differ min1: " << min1 << "!=" << idx1 << ". min2: " << min2 << "!=" << idx2 << ". Value: " << min << "!=" << value << endl;
      exit(1);
    }
    //return mytree->serialize_tree(activeRows[index1],activeRows[index2], distance);  
  } else {
    findMin();    
    //    cout << "TEST:    JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << min << " clusterCount=" << clusterCount << "\n";
    if(idx1 == min1 && idx2 == min2){
      if((value - min) * (value - min) > 0.00001){
	printf("%.20f  %.20f \n",value,min);
	cout << "ERROR: indexes match but value doesn't: " << value << "!=" << min << endl;  
	exit(1);
      } 
    } else if((value - min) * (value - min) > 0.00001){
      printf("%.20f  %.20f \n",value,min);
      cout << "ERROR: both value and indexes differ min1: " << min1 << "!=" << idx1 << ". min2: " << min2 << "!=" << idx2 << ". Value: " << min << "!=" << value << endl;
      cout << matrix[970][943] << "-" << separations[970] << "-" << separations[943] << " = " << matrix[970][943] - separations[970] - separations[943] << endl; 
      exit(1);
    } else {
      min1 = idx1;
      min2 = idx2;
    }  
  }
  mergeMinNodes();
  updateMatrix();
}

void testNJ::stepNoVal(int idx1, int idx2){
  //  cout << idx2 << " " << idx1 << endl;
  if(clusterCount < 4){
    return;
  }
  if( clusterCount == 2){
    // finish by joining the two remaining clusters
    min1 = -1;
    min2 = -1;
    // find the last nodes
    for(int i = 0; i < matrixSize; i++){
      if(activeRows[i] != -1){
	if(min1 == -1){
	  min1 = i;
	} else {
	  min2 = i;	
	  break;
	}            
      }
    }  
    //    double distance = matrix[min1][min2];
    if(!((idx1 == min1 && idx2 == min2) || (idx2 == min1 && idx1 == min2))){   
      double checkVal = matrix[idx1][idx2] - separations[idx1] - separations[idx2];
      if(checkVal > min){
	cout << "ERROR: bad join min1: " << min1 << "!=" << idx1 << ". min2: " << min2 << "!=" << idx2 << ". Value: " << min << "!=" << checkVal << endl << endl;
      }      
      exit(1);
    }
    //return mytree->serialize_tree(activeRows[index1],activeRows[index2], distance);  
  } else {
    findMin();    
    //cout << "TEST:    JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << min << " clusterCount=" << clusterCount << "\n";
    if(!((idx1 == min1 && idx2 == min2) || (idx2 == min1 && idx1 == min2))){   
      //check if this join is ok
      double checkVal = matrix[idx1][idx2] - separations[idx1] - separations[idx2];
      if((checkVal - min) * (checkVal - min) > 0.0001){
	printf("%.20f  %.20f \n",checkVal,min);
	cout << "ERROR: bad join min1: " << min1 << "!=" << idx1 << ". min2: " << min2 << "!=" << idx2 << ". Value: " << min << "!=" << checkVal << endl << endl;
      }
      exit(1);
    } else {
      min1 = idx1;
      min2 = idx2;
    }
  }
  mergeMinNodes();
  updateMatrix();
}


void testNJ::updateMatrix(){    
  distType newSeparationsum = 0;
  distType mutualDistance = matrix[min1][min2];
  distType* row1 = matrix[min1];
  distType* row2 = matrix[min2];
  for(int i = 0; i < matrixSize; i++){
    if(i == min1 || i == min2 || activeRows[i] == -1){
      row1[i] = 0;
    } else {
      distType val1 = row1[i];
      distType val2 = row2[i];
      distType dist = (val1 + val2 - mutualDistance) / ((distType)2.0);
      newSeparationsum += dist;
      // update the separationsum of cluster i.
      separationsums[i] += (dist - val1 - val2);
      
      separations[i] = separationsums[i] / (clusterCount -2); 
      row1[i] = dist;
      matrix[i][min1] = dist;
    }
  }
  separationsums[min1] = newSeparationsum;
  separations[min1] = newSeparationsum / (clusterCount - 2);
  separationsums[min2] = 0;
  activeRows[min2] = -1;
  activeRows[min1] =  nextId++ ;  
}

void testNJ::mergeMinNodes(){
  // calculate distances
  double dist = matrix[min1][min2];
  double sep1 = separations[min1];
  double sep2 = separations[min2];
  double dist1 = (0.5 * dist) + (0.5 * (sep1 - sep2));
  double dist2 = (0.5 * dist) + (0.5 * (sep2 - sep1));
  // update tree
  mytree->addInternalNode(dist1, dist2, activeRows[min1], activeRows[min2]);
  clusterCount--;
}

void testNJ::initialize(){
  mytree = new polytree(matrixSize, reader->getSequenceNames());  
  //calculate initial seperation rows
  for(int i = 0; i < matrixSize; i++){
    distType sum = 0;
    for(int j = 0; j < matrixSize; j++){
      sum += matrix[i][j];
    }
    separationsums[i] = sum;
    separations[i] = sum / (clusterCount - 2); 
    activeRows[i] = i;
  }
}

void testNJ::findMin() {
  min1 = -1;
  min2 = -1;
  min = DBL_MAX;
  for (int i = 0; i < matrixSize; i++) {    
    if(activeRows[i] != -1){
      distType* row = matrix[i];
      double sep1 = separations[i];
      for(int j = 0; j < matrixSize; j++){
	if(activeRows[j] != -1 && i != j){
	  double sep2 = separations[j];
	  double val = row[j] - sep1 - sep2;
	  //	    cout << "TEST:  [" << i << "," << j << "] " << row[j] << "-" << sep1 << "-" << sep2 << "=" << val << "    min=" << min << endl;
	  if(val < min){
	    // new minimum
	    min1 = i;
	    min2 = j;
	    min = val;
	  }					
	}
      }
    }
  }
  //  cout << "TEST:    JOINING min1=" << min1 << " min2=" << min2 << " global_min=" << min << " clusterCount=" << clusterCount << "\n";
}

testNJ::~testNJ(void)
{
}
