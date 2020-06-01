#ifndef CLUSTER_PAIR_H
#define CLUSTER_PAIR_H

struct cluster_pair{  
public:
  unsigned int id; /* index of the other cluster*/
  distType distance; /* distance between clusters */
  bool operator<(const cluster_pair& a) const{
    return distance < a.distance;
  }  
};  

#endif
