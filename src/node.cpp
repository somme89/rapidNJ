#include "stdinclude.h"
#include "node.h"
using namespace std;

/* This class represents a node in the final unrooted tree.
 * It has the nessasary logic to serialize the tree to a Newick formatted string representation. 
*/

/* Creates a internal node*/
node::node() {	
  edgeDist = new distType[3];
  edges = new node*[3];
  edgeCount = 0;
}

/* creates a named leaf node */
node::node(string name) {
  edgeDist = new distType[1];
  edges = new node*[1];
  node::name = name;
  edgeCount = 0;
}

/* Connect a node to this node */
void node::addEdge(node* n, distType distance){
  if(n == this){
    cout << "cannot add node to self \n";
    exit(0);
  }
  edgeDist[edgeCount] = distance;
  
  edges[edgeCount] = n;
  edgeCount++;
}

/* serialize the three to newick format*/
string node::serializeTree(){
  if(edgeCount == 3){
    //internal node
    return serializeNode(NULL) + ";";
  } else {
    if(edgeCount == 0){
      // only one node, this one.
      return name + ";";
    } else if(edgeCount == 1 && edges[0]->edgeCount == 1){
      // only two nodes.
      return "(" + name + "," + edges[0]->name + ");";
    } else {
      // leaf, use the connected internal node to start the recursion
      return edges[0]->serializeTree();
    }
  }
  return "ERROR - could not serialize the tree";
}

/* Recursive helper method for serializing the three*/
string node::serializeNode(node* n){
  if(edgeCount == 2){
    cout << "ERROR - only 2 edges found\n";
  }
  char buffer [50];
  if(edgeCount == 3){
    // internal node
    string s = "(";
    for(int i = 0; i < edgeCount; i++){
      node* edge = edges[i];
      if(edge != n){
	s += edge->serializeNode(this);
	s += ":";				
	int length = sprintf(buffer,"%g", edgeDist[i]);
	s.append(buffer,length);
	s += ",";
      }
    }
    s.replace(s.length()-1,1,")");
    return s;
  } else {
    return "'" + name + "'";
  }
  
}

node::~node(void)
{
}
