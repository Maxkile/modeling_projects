#include "vmo.hpp"

// Compute eucledian norm
double vmo::norm(const std::vector<double> &vec, NodesInfo *nodesinfo) { return sqrt(dot(vec, vec, nodesinfo)); }
