#ifndef NODE_H
#define NODE_H
#include <string>

class node
{
public:
  node();
  node(std::string name);
  void addEdge(node* n, distType distance);
  ~node(void);
  
  std::string name;
  distType* edgeDist;
  node** edges;
  std::string serializeTree();
  std::string serializeNode(node* n);
  int edgeCount;

private:
    
};

#endif
