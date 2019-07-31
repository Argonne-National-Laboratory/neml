#ifndef SPHERE_H
#define SPHERE_H

/*
 *  The associated .cxx file contains the efficient spherical designs,
 *  in spherical coordinates, 
 *  from http://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/index.html
 */

#include <vector>

#define MAX_S2_DEGREE 180
#define MAX_SO3_DEGREE 16 // 31

namespace neml {

extern const std::vector<std::vector<double>> S2_theta;
extern const std::vector<std::vector<double>> S2_phi;
extern const std::vector<std::vector<std::vector<double>>> SO3_pts;

}

#endif
