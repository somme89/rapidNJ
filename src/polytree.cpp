/*This class is a very efficient implementation af a bifurcating phylogenetic tree.
* It can be serialised to Newick format. 
*/

#include "stdinclude.h"
#include "polytree.h"

using namespace std;

/* 
* size = number of leafs
*/
polytree::polytree(int size, vector<string>* sequenceNames){
  distances = new double[size*2];
  left_indexes = new int[size];
  right_indexes = new int[size];
  leaf_names = new string[size];
  polytree::size = size;
  parent_indices = new int[size*2];
  bootstrap_counts = new int[size];
  bootstrap_replicate_count = 0;
  test = false;
  hmap = NULL;
  for (int i = 0; i < size*2; i++) {
    parent_indices[i] = -1;
    if(i < size) {
      bootstrap_counts[i] = 0;
    }
  }
  current_index = 0;
  leaf_index = 0;  
  if(sequenceNames != NULL) {
    for(int k = 0; k < size; k++){
      addLeaf(sequenceNames->at(k));
    }
  }
}

void polytree::addLeaf(string name){
  leaf_names[leaf_index] = name;
  leaf_index++;
}

void polytree::addInternalNode(double left_dist, double right_dist, int left_index, int right_index){  
  left_indexes[current_index] = left_index;
  right_indexes[current_index] = right_index;    
  distances[left_index] = left_dist;
  distances[right_index] = right_dist;
  parent_indices[left_index] = current_index + size;
  parent_indices[right_index] = current_index + size;
  current_index++;
}

int polytree::getParentIndex(int child_index){
  return parent_indices[child_index];
}

int polytree::getLeftChildIndex(int parentIndex) {
  return left_indexes[parentIndex - size];
}

int polytree::getRightChildIndex(int parentIndex) {  
  return right_indexes[parentIndex - size];
}

int polytree::getSibling(int index) {
  int parent = getParentIndex(index);
  if(parent == -1){
    if(index == s_index_left && s_index_right >= size) {
      parent = s_index_right;
    } else if (s_index_left >= size){
      parent = s_index_left;
    } else {
      return -1;
    }
  }
  if(getLeftChildIndex(parent) == index) {
    return getRightChildIndex(parent);
  } else {
    return getLeftChildIndex(parent);
  }
}

void polytree::InitLeafLists(vector<int>** list) {
  for(unsigned int i = 0; i < size*2; i++) {
    list[i] = NULL;
  }
  for(unsigned int i = 0; i < size; i++) {    
    list[i] = new vector<int>();    
    list[i]->push_back(i);
  }
}

/** Combine the two lists and ensure that the new list is sorted*/
vector<int>* combineLeafLists(vector<int>* l1, vector<int>* l2) {
  vector<int>* retVal = new vector<int>;
  unsigned int idx1 = 0;
  unsigned int idx2 = 0;
  while(retVal->size() < l1->size() + l2->size()) {
    if(idx1 == l1->size()) {
      retVal->push_back(l2->at(idx2));
      idx2++;
    } else if(idx2 == l2->size()) {
      retVal->push_back(l1->at(idx1));
      idx1++;
    }
    else if(l1->at(idx1) < l2->at(idx2)) {
      retVal->push_back(l1->at(idx1));
      idx1++;
    } else {
      retVal->push_back(l2->at(idx2));
      idx2++;
    }
  }  
  return retVal;
}

unsigned int polytree::increaseLeafListsSize(vector<int>** leafLists, unsigned int newMaxSize) {
  unsigned int mergeCount = 0;
  for(unsigned int i = 0; i < size*2; i++) {
    if(leafLists[i] == NULL) {
      continue;
    }
    int parent = getParentIndex(i);
    if(parent == -1) {
      continue;
    }
    int sibling = getSibling(i);
    if(leafLists[sibling] != NULL) {
      if(leafLists[i]->size() + leafLists[sibling]->size() <= newMaxSize) {
        //Combine them
        vector<int>* newList = combineLeafLists(leafLists[i], leafLists[sibling]);
        leafLists[parent] = newList;
        delete leafLists[i];
        leafLists[i] = NULL;
        delete leafLists[sibling];
        leafLists[sibling] = NULL;
        mergeCount++;
      }
    }
  }
  return mergeCount;
}

int isEqualLists(vector<int>* l1, vector<int>* l2) {
  if(l1 == NULL || l2 == NULL) {
    return 0;
  }
  if(l1->size() != l2->size()) {
    return 0;
  }
  for(unsigned int i = 0; i < l1->size(); i++) {
    if(l1->at(i) != l2->at(i)) {
      return 0;
    } 
  }
  return 1;
}

int isEqualLists2(vector<int>* l1, vector<int>* l2) {
  if(l1 == NULL || l2 == NULL) {
    return 0;
  }
  if(l1->size() != l2->size()) {
    return 0;
  }
  for(unsigned int i = 0; i < l1->size(); i++) {
    if(l1->at(i) != l2->at(i)) {
      return 0;
    } 
  }
  return 1;
}

void polytree::setlastParentIndicies(){
  if(s_index_right >= size) {
    parent_indices[s_index_left] = s_index_right;
  }
  if(s_index_left >= size) {
    parent_indices[s_index_right] = s_index_left;
  }
}

void polytree::getNextElements(stack<int>* nextElements, int* nodeUsed, stack<int>* tempElements, vector<int>* section){
  while(!tempElements->empty()) {
    int node = tempElements->top();
    tempElements->pop();
    int p = getParentIndex(node);
    if(p >= size) {
      //Internal node
      if(nodeUsed[p] == 0) {
        nextElements->push(p);
      }
    } else{
      //Leaf node
      section->push_back(p);
    }
    nodeUsed[p] = 1;
    int c1 = getRightChildIndex(node);
    if(c1 >= size) {
      if(nodeUsed[c1] == 0) {
        nextElements->push(c1);
      }
    } else{
      section->push_back(c1);
    }
    nodeUsed[c1] = 1;
    int c2 = getLeftChildIndex(node);
    if(c2 >= size) {
      if(nodeUsed[c2] == 0) {
        nextElements->push(c2);
      }
    } else{
      section->push_back(c2);
    }
    nodeUsed[c2] = 1;
  }
}

bisection* polytree::getBisectionsForNodePair(int n1, int n2) {
  int* nodeUsed = new int[size*2];
  stack<int>* nextElements = new stack<int>();
  stack<int>* tempElements = new stack<int>();
  for(int i = 0; i < size*2; i++){
    nodeUsed[i] = 0;
  }
  nodeUsed[n1] = 1;
  nodeUsed[n2] = 1;
  bisection* bs = new bisection(n1, n2);

  //traverse the tree and discover members of each bipartition
  tempElements->push(n1);
  while(!tempElements->empty()) {
    getNextElements(nextElements, nodeUsed, tempElements, bs->section1);
    stack<int>* t = nextElements;
    nextElements = tempElements;
    tempElements = t;
  }
  tempElements->push(n2);
  while(!tempElements->empty()) {
    getNextElements(nextElements, nodeUsed, tempElements, bs->section2);
    stack<int>* t = nextElements;
    nextElements = tempElements;
    tempElements = t;
  }
  delete tempElements;
  delete nextElements;
  delete[] nodeUsed;
  return bs;
}

void polytree::findBisections(int i, int candidate, int* nodeUseCount, map<int, vector<bisection*>*>* hmap, bool useIndexAsHash){
  bisection* bisec = getBisectionsForNodePair(i, candidate);
  int h = bisec->hash();
  if(useIndexAsHash){
    h = candidate;
  }
  if(hmap->find(h) == hmap->end()) {
    vector<bisection*>* v = new vector<bisection*>();
    v->push_back(bisec);
    hmap->insert(make_pair(h,v));
  } else {
    hmap->at(h)->push_back(bisec);
  }
  nodeUseCount[candidate]++; 
}

map<int, vector<bisection*>*>* polytree::findAllBisections(bool useIndexAsHash){
  int* nodeUseCount = new int[size*2];
  map<int, vector<bisection*>*>* hmap = new map<int, vector<bisection*>*>();
  for(int i = size; i < size*2; i++){
    nodeUseCount[i] = 0;
  }
  //Do a breath first search to enable indexing of the initial tree using inner nodes.
  vector<unsigned int>* currentNodes = new vector<unsigned int>;
  vector<unsigned int>* nextNodes = new vector<unsigned int>;

  currentNodes->push_back(size*2-3); //start at last inner node
  while(currentNodes->size() != 0){
    unsigned int n1 = currentNodes->back();
    currentNodes->pop_back();
    if(n1 < size) {
      cerr << "ERROR: Node with id " << n1 << " is a leaf node" << endl;
      exit(1);
    }
    nodeUseCount[n1] = 3;
    int n2 = getRightChildIndex(n1);
    if(n2 >= size && nodeUseCount[n2] < 3) {
      findBisections(n1, n2, nodeUseCount, hmap, useIndexAsHash);
      nextNodes->push_back(n2);
    }
    n2 = getLeftChildIndex(n1);
    if(n2 >= size && nodeUseCount[n2] < 3) {
      findBisections(n1, n2, nodeUseCount, hmap, useIndexAsHash);
      nextNodes->push_back(n2);
    }
    n2 = getParentIndex(n1);
    if(n2 >= size && nodeUseCount[n2] < 3) {
      findBisections(n1, n2, nodeUseCount, hmap, useIndexAsHash);
      nextNodes->push_back(n2);
    }
    if(currentNodes->size() == 0) {
      vector<unsigned int>* temp = currentNodes;
      currentNodes = nextNodes;
      nextNodes = temp;
    }
  }
  delete[] nodeUseCount;
  delete currentNodes;
  delete nextNodes;
  return hmap;
}

void polytree::updateBSValues(unsigned int n1, unsigned int n2, map<int, vector<bisection*>*>* hmap, map<int, vector<bisection*>*>* hmapOther){
  map<int, vector<bisection*>*>::iterator it = hmap->find(n2);
  if(it == hmap->end()) {
    return;
  }
  vector<bisection*>* n2_bisecs = it->second;
  if(n2_bisecs->size() > 1) {
    cout << "ERROR: multiple bi-sections assigned to internal node with id: " << n2 << endl;
    exit(1);
  }
  bisection* n2_bisec = n2_bisecs->at(0);
  map<int, vector<bisection*>*>::iterator itOther = hmapOther->find(n2_bisecs->back()->hash());
  if(itOther == hmapOther->end()) {
    return;
  }
  vector<bisection*>* bisecsOther = itOther->second;
  for(int i = 0; i < bisecsOther->size(); i++) {
    if(n2_bisec->equals(bisecsOther->at(i))){
      if(bootstrap_counts[n2-size] == -1) {
        bootstrap_counts[n2-size] = 0;
      }
      bootstrap_counts[n2-size]++;
    }
  }  
}

void polytree::compareTreeBootstrap(polytree* tree) {
  bootstrap_replicate_count++;
  if(tree->size != size) {
    cerr << "Cannot compare trees with different sizes" << endl;
    exit(1);
  }
  //Find and store all bisections in the two trees 
  if(hmap == NULL) {
    hmap = findAllBisections(true);
  }
  map<int, vector<bisection*>*>* hmapOther = tree->findAllBisections(false);
  //Do a breath first traversel of the original tree and search for each bisection in the original tree in the new tree.
  vector<unsigned int>* currentNodes = new vector<unsigned int>;
  vector<unsigned int>* nextNodes = new vector<unsigned int>;

  int* visited = new int[size*2];
  for(int i = 0; i < size*2; i++){
    visited[i] = 0;
  }
  currentNodes->push_back(size*2-3); //start at last inner node
  while(currentNodes->size() != 0){
    unsigned int n1 = currentNodes->back();
    currentNodes->pop_back();
    if(n1 < size) {
      cerr << "ERROR: Node with id " << n1 << " is a leaf node" << endl;
      exit(1);
    }
    visited[n1] = 1;
    int n2 = getRightChildIndex(n1);
    if(n2 >= size && visited[n2] == 0) {
      nextNodes->push_back(n2);
      updateBSValues(n1, n2, hmap, hmapOther);
    }
    n2 = getLeftChildIndex(n1);
    if(n2 >= size && visited[n2] == 0) {
      nextNodes->push_back(n2);
      updateBSValues(n1, n2, hmap, hmapOther);
    }
    n2 = getParentIndex(n1);
    if(n2 >= size && visited[n2] == 0) {
      nextNodes->push_back(n2);
      updateBSValues(n1, n2, hmap, hmapOther);
    }
    if(currentNodes->size() == 0) {
      vector<unsigned int>* temp = currentNodes;
      currentNodes = nextNodes;
      nextNodes = temp;
    }
  }
  delete[] visited;
  delete currentNodes;
  delete nextNodes;
  for ( map<int, vector<bisection*>*>::iterator it = hmapOther->begin(); it != hmapOther->end(); ++it ){
    vector<bisection*>* v = it->second;
    for( vector<bisection*>::iterator it2 = v->begin(); it2 != v->end(); ++it2){
      delete *it2;
    }
    delete v;
  } 
  delete hmapOther;
}

void polytree::serialize_tree(ostream &out){  
  /**for (int i = 0; i < size; i++) {
      cout << i << " " << bootstrap_counts[i] << endl;
    }
  */
  
  out << "(";
  if(s_index_left < size){
    out << "'" << leaf_names[s_index_left] << "'";
    out << ":";
    int length = sprintf(buffer,"%.5g", s_dist);
    string s = "";
    s.append(buffer,length);
    out << s;
  } else {
    if(s_index_right >= size){
      out << "(";
    }
    serialize_node(out, left_indexes[s_index_left - size] , right_indexes[s_index_left - size], s_index_left);
    if(s_index_right >= size){
      out << ")";
      if(bootstrap_replicate_count > 0 && s_index_left < size*2-3) {        
        int bs = (int) (((double)bootstrap_counts[s_index_left-size] / (double) bootstrap_replicate_count) * 100);
        out << bs;
      }
      out << ":";
      int length = sprintf(buffer,"%.5g", s_dist);
      string s = "";
      s.append(buffer,length);
      out << s;
    }
  }
  out << ",";
  if(s_index_right < size){
    out << "'" << leaf_names[s_index_right] << "'";
    out << ":";
    int length = sprintf(buffer,"%.5g", s_dist);
    string s = "";
    s.append(buffer,length);
    out << s;
  } else { 
    serialize_node(out, left_indexes[s_index_right - size], right_indexes[s_index_right - size],s_index_right);    
  }
  out << ")";

  if(bootstrap_replicate_count > 0) {
    if (s_index_right >= size && s_index_right < size*2-3) {
      int bs = (int) (((double)bootstrap_counts[s_index_right-size] / (double) bootstrap_replicate_count) * 100);
      out << bs;
    } else if (s_index_left >= size && s_index_left < size*2-3){
      int bs = (int) (((double)bootstrap_counts[s_index_left-size] / (double) bootstrap_replicate_count) * 100);
      out << bs;
    }    
  }
  out << ";" << endl;
}

void polytree::serialize_node(ostream &out, unsigned int left_index, unsigned int right_index, int index){
  if(left_index < size){
    // leaf node
    out << "'" << leaf_names[left_index] << "'";
    out << ":";
    int length = sprintf(buffer,"%.5g", distances[left_index]);
    string s = "";
    s.append(buffer,length);
    out << s;
  } else {
    // serialize the left tree recursively
    out << "(";
    serialize_node(out, left_indexes[left_index - size], right_indexes[left_index-size], left_index);
    out << ")";
    if(bootstrap_replicate_count > 0 && left_index < size*2-3) {
      int bs = (int) (((double)bootstrap_counts[left_index-size] / (double) bootstrap_replicate_count) * 100);
      out << bs;
    }
    out << ":";
    int length = sprintf(buffer,"%.5g", distances[left_index]);
    string s = "";
    s.append(buffer,length);
    out << s;
  }
  out << ",";
  if(right_index < size){
    // leaf node
    out << "'" << leaf_names[right_index] << "'";
    out << ":";
    int length = sprintf(buffer,"%.5g", distances[right_index]);
    string s = "";
    s.append(buffer,length);
    out << s;
  }
  else{
    // serialize the right tree recursively
    out << "(";
    serialize_node(out, left_indexes[right_index - size], right_indexes[right_index-size], right_index);
    out << ")";
    if(bootstrap_replicate_count > 0 && right_index < size*2-3) {
      int bs = (int) (((double)bootstrap_counts[right_index-size] / (double) bootstrap_replicate_count) * 100);
      out << bs;      
    }
    out << ":";
    int length = sprintf(buffer,"%.5g", distances[right_index]);
    string s = "";
    s.append(buffer,length);
    out << s;
  }
}

/*Set the indices of the last two clusters that have are joined.*/
void polytree::set_serialization_indices(int left, int right, distType dist) {
  s_index_left = left;
  s_index_right = right;
  s_dist = dist;
  //Parent indicies have to be set to allow traversal of the tree in bootstrapping.
  parent_indices[s_index_left] = s_index_right;
  parent_indices[s_index_right] = s_index_left;

}

polytree::~polytree(void){
  delete[] distances;
  delete[] left_indexes;
  delete[] right_indexes;
  delete[] leaf_names;
  delete[] parent_indices;
  delete[] bootstrap_counts;
  if(hmap != NULL){
    delete hmap;
  }
}
