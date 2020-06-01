#include "DistanceEstimate.hpp"
#include "bitDistanceGap.hpp"
#include "bitDistanceProtein.hpp"
#include "dataloader.hpp"
#include "simpleDistanceCalculator.hpp"

DistanceEstimate::DistanceEstimate(dataloader* loader) {
  DistanceEstimate::type = loader->type;
  DistanceEstimate::fastDist = loader->fastdist;
  DistanceEstimate::loader = loader;
}

DistanceEstimate::~DistanceEstimate(void) {

}
