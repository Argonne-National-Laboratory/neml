#include "math/projections.h"

namespace neml {

RankFour normal_projection(const Vector & n)
{
  RankFour S;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        for (size_t l = 0; l < 3; l++) {
          S(i,j,k,l) += n(k) * n(l) * n(i) * n(j);
        }
      }
    }
  }

  return S;
}

SymSymR4 normal_projection_ss(const Vector & n)
{
  double sq2 = std::sqrt(2.0);

  SymSymR4 SS;

	SS(0,0) = n(0)*n(0)*n(0)*n(0);
	SS(0,1) = n(0)*n(0)*n(1)*n(1);
	SS(0,2) = n(0)*n(0)*n(2)*n(2);
	SS(0,3) = sq2*n(0)*n(0)*n(1)*n(2);
	SS(0,4) = sq2*n(0)*n(0)*n(0)*n(2);
	SS(0,5) = sq2*n(0)*n(0)*n(0)*n(1);
	SS(1,0) = n(0)*n(0)*n(1)*n(1);
	SS(1,1) = n(1)*n(1)*n(1)*n(1);
	SS(1,2) = n(1)*n(1)*n(2)*n(2);
	SS(1,3) = sq2*n(1)*n(1)*n(1)*n(2);
	SS(1,4) = sq2*n(0)*n(1)*n(1)*n(2);
	SS(1,5) = sq2*n(0)*n(1)*n(1)*n(1);
	SS(2,0) = n(0)*n(0)*n(2)*n(2);
	SS(2,1) = n(1)*n(1)*n(2)*n(2);
	SS(2,2) = n(2)*n(2)*n(2)*n(2);
	SS(2,3) = sq2*n(1)*n(2)*n(2)*n(2);
	SS(2,4) = sq2*n(0)*n(2)*n(2)*n(2);
	SS(2,5) = sq2*n(0)*n(1)*n(2)*n(2);
	SS(3,0) = sq2*n(0)*n(0)*n(1)*n(2);
	SS(3,1) = sq2*n(1)*n(1)*n(1)*n(2);
	SS(3,2) = sq2*n(1)*n(2)*n(2)*n(2);
	SS(3,3) = 2*n(1)*n(1)*n(2)*n(2);
	SS(3,4) = 2*n(0)*n(1)*n(2)*n(2);
	SS(3,5) = 2*n(0)*n(1)*n(1)*n(2);
	SS(4,0) = sq2*n(0)*n(0)*n(0)*n(2);
	SS(4,1) = sq2*n(0)*n(1)*n(1)*n(2);
	SS(4,2) = sq2*n(0)*n(2)*n(2)*n(2);
	SS(4,3) = 2*n(0)*n(1)*n(2)*n(2);
	SS(4,4) = 2*n(0)*n(0)*n(2)*n(2);
	SS(4,5) = 2*n(0)*n(0)*n(1)*n(2);
	SS(5,0) = sq2*n(0)*n(0)*n(0)*n(1);
	SS(5,1) = sq2*n(0)*n(1)*n(1)*n(1);
	SS(5,2) = sq2*n(0)*n(1)*n(2)*n(2);
	SS(5,3) = 2*n(0)*n(1)*n(1)*n(2);
	SS(5,4) = 2*n(0)*n(0)*n(1)*n(2);
	SS(5,5) = 2*n(0)*n(0)*n(1)*n(1);

  return SS;
}

RankFour shear_projection(const Vector & n)
{
  RankFour S;
  RankTwo I = RankTwo::id();

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        for (size_t l = 0; l < 3; l++) {
          S(i,j,k,l) = (I(i,k) - n(i) * n(k)) * n(j) * n(l);
        }
      }
    }
  }
 
  return S;
}

SymSymR4 shear_projection_ss(const Vector & n)
{
  SymSymR4 SS;

  double sq2 = std::sqrt(2.0);

	SS(0,0) = 2*n(0)*n(0)*(1 - n(0)*n(0));
	SS(0,1) = -2*n(0)*n(0)*n(1)*n(1);
	SS(0,2) = -2*n(0)*n(0)*n(2)*n(2);
	SS(0,3) = -2*sq2*n(0)*n(0)*n(1)*n(2);
	SS(0,4) = sq2*n(0)*n(2)*(1 - 2*n(0)*n(0));
	SS(0,5) = sq2*n(0)*n(1)*(1 - 2*n(0)*n(0));
	SS(1,0) = -2*n(0)*n(0)*n(1)*n(1);
	SS(1,1) = 2*n(1)*n(1)*(1 - n(1)*n(1));
	SS(1,2) = -2*n(1)*n(1)*n(2)*n(2);
	SS(1,3) = sq2*n(1)*n(2)*(1 - 2*n(1)*n(1));
	SS(1,4) = -2*sq2*n(0)*n(1)*n(1)*n(2);
	SS(1,5) = sq2*n(0)*n(1)*(1 - 2*n(1)*n(1));
	SS(2,0) = -2*n(0)*n(0)*n(2)*n(2);
	SS(2,1) = -2*n(1)*n(1)*n(2)*n(2);
	SS(2,2) = 2*n(2)*n(2)*(1 - n(2)*n(2));
	SS(2,3) = sq2*n(1)*n(2)*(1 - 2*n(2)*n(2));
	SS(2,4) = sq2*n(0)*n(2)*(1 - 2*n(2)*n(2));
	SS(2,5) = -2*sq2*n(0)*n(1)*n(2)*n(2);
	SS(3,0) = -2*sq2*n(0)*n(0)*n(1)*n(2);
	SS(3,1) = sq2*n(1)*n(2)*(1 - 2*n(1)*n(1));
	SS(3,2) = sq2*n(1)*n(2)*(1 - 2*n(2)*n(2));
	SS(3,3) = -4*n(1)*n(1)*n(2)*n(2) + n(1)*n(1) + n(2)*n(2);
	SS(3,4) = n(0)*n(1)*(1 - 4*n(2)*n(2));
	SS(3,5) = n(0)*n(2)*(1 - 4*n(1)*n(1));
	SS(4,0) = sq2*n(0)*n(2)*(1 - 2*n(0)*n(0));
	SS(4,1) = -2*sq2*n(0)*n(1)*n(1)*n(2);
	SS(4,2) = sq2*n(0)*n(2)*(1 - 2*n(2)*n(2));
	SS(4,3) = n(0)*n(1)*(1 - 4*n(2)*n(2));
	SS(4,4) = -4*n(0)*n(0)*n(2)*n(2) + n(0)*n(0) + n(2)*n(2);
	SS(4,5) = n(1)*n(2)*(1 - 4*n(0)*n(0));
	SS(5,0) = sq2*n(0)*n(1)*(1 - 2*n(0)*n(0));
	SS(5,1) = sq2*n(0)*n(1)*(1 - 2*n(1)*n(1));
	SS(5,2) = -2*sq2*n(0)*n(1)*n(2)*n(2);
	SS(5,3) = n(0)*n(2)*(1 - 4*n(1)*n(1));
	SS(5,4) = n(1)*n(2)*(1 - 4*n(0)*n(0));
	SS(5,5) = -4*n(0)*n(0)*n(1)*n(1) + n(0)*n(0) + n(1)*n(1);

  return SS;
}

} // namespace neml
