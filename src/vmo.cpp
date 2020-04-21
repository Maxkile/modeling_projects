#include "vmo.hpp"

// Compute eucledian norm
double vmo::norm(const std::vector<double>& vec){
    return sqrt(dot(vec,vec));
}
