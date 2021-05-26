#pragma once

#include "tensors.h"

namespace neml {

/// Projection operator onto the plane and normal direction
RankFour normal_projection(const Vector & n);
/// Plane and normal direction with minor symmetry on both sides
SymSymR4 normal_projection_ss(const Vector & n);

/// Projection operator onto the shear plane
RankFour shear_projection(const Vector & n);
/// Shear projection with minor symmetry on both sides
SymSymR4 shear_projection_ss(const Vector & n);


} // namespace neml
