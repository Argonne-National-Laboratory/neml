#include "nemlmath.h"

#include "../nemlerror.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace neml {

int transform_fourth(const double * const D, const double * const W, double * const M)
{
	M[0] = D[0];
	M[1] = sqrt(2)*D[5]/2 - W[2]/2;
	M[2] = sqrt(2)*D[4]/2 + W[1]/2;
	M[3] = sqrt(2)*D[5]/2 + W[2]/2;
	M[4] = D[1];
	M[5] = sqrt(2)*D[3]/2 - W[0]/2;
	M[6] = sqrt(2)*D[4]/2 - W[1]/2;
	M[7] = sqrt(2)*D[3]/2 + W[0]/2;
	M[8] = D[2];
	M[9] = sqrt(2)*D[30]/2;
	M[10] = D[35]/2 - sqrt(2)*W[17]/4;
	M[11] = D[34]/2 + sqrt(2)*W[16]/4;
	M[12] = D[35]/2 + sqrt(2)*W[17]/4;
	M[13] = sqrt(2)*D[31]/2;
	M[14] = D[33]/2 - sqrt(2)*W[15]/4;
	M[15] = D[34]/2 - sqrt(2)*W[16]/4;
	M[16] = D[33]/2 + sqrt(2)*W[15]/4;
	M[17] = sqrt(2)*D[32]/2;
	M[18] = sqrt(2)*D[24]/2;
	M[19] = D[29]/2 - sqrt(2)*W[14]/4;
	M[20] = D[28]/2 + sqrt(2)*W[13]/4;
	M[21] = D[29]/2 + sqrt(2)*W[14]/4;
	M[22] = sqrt(2)*D[25]/2;
	M[23] = D[27]/2 - sqrt(2)*W[12]/4;
	M[24] = D[28]/2 - sqrt(2)*W[13]/4;
	M[25] = D[27]/2 + sqrt(2)*W[12]/4;
	M[26] = sqrt(2)*D[26]/2;
	M[27] = sqrt(2)*D[30]/2;
	M[28] = D[35]/2 - sqrt(2)*W[17]/4;
	M[29] = D[34]/2 + sqrt(2)*W[16]/4;
	M[30] = D[35]/2 + sqrt(2)*W[17]/4;
	M[31] = sqrt(2)*D[31]/2;
	M[32] = D[33]/2 - sqrt(2)*W[15]/4;
	M[33] = D[34]/2 - sqrt(2)*W[16]/4;
	M[34] = D[33]/2 + sqrt(2)*W[15]/4;
	M[35] = sqrt(2)*D[32]/2;
	M[36] = D[6];
	M[37] = sqrt(2)*D[11]/2 - W[5]/2;
	M[38] = sqrt(2)*D[10]/2 + W[4]/2;
	M[39] = sqrt(2)*D[11]/2 + W[5]/2;
	M[40] = D[7];
	M[41] = sqrt(2)*D[9]/2 - W[3]/2;
	M[42] = sqrt(2)*D[10]/2 - W[4]/2;
	M[43] = sqrt(2)*D[9]/2 + W[3]/2;
	M[44] = D[8];
	M[45] = sqrt(2)*D[18]/2;
	M[46] = D[23]/2 - sqrt(2)*W[11]/4;
	M[47] = D[22]/2 + sqrt(2)*W[10]/4;
	M[48] = D[23]/2 + sqrt(2)*W[11]/4;
	M[49] = sqrt(2)*D[19]/2;
	M[50] = D[21]/2 - sqrt(2)*W[9]/4;
	M[51] = D[22]/2 - sqrt(2)*W[10]/4;
	M[52] = D[21]/2 + sqrt(2)*W[9]/4;
	M[53] = sqrt(2)*D[20]/2;
	M[54] = sqrt(2)*D[24]/2;
	M[55] = D[29]/2 - sqrt(2)*W[14]/4;
	M[56] = D[28]/2 + sqrt(2)*W[13]/4;
	M[57] = D[29]/2 + sqrt(2)*W[14]/4;
	M[58] = sqrt(2)*D[25]/2;
	M[59] = D[27]/2 - sqrt(2)*W[12]/4;
	M[60] = D[28]/2 - sqrt(2)*W[13]/4;
	M[61] = D[27]/2 + sqrt(2)*W[12]/4;
	M[62] = sqrt(2)*D[26]/2;
	M[63] = sqrt(2)*D[18]/2;
	M[64] = D[23]/2 - sqrt(2)*W[11]/4;
	M[65] = D[22]/2 + sqrt(2)*W[10]/4;
	M[66] = D[23]/2 + sqrt(2)*W[11]/4;
	M[67] = sqrt(2)*D[19]/2;
	M[68] = D[21]/2 - sqrt(2)*W[9]/4;
	M[69] = D[22]/2 - sqrt(2)*W[10]/4;
	M[70] = D[21]/2 + sqrt(2)*W[9]/4;
	M[71] = sqrt(2)*D[20]/2;
	M[72] = D[12];
	M[73] = sqrt(2)*D[17]/2 - W[8]/2;
	M[74] = sqrt(2)*D[16]/2 + W[7]/2;
	M[75] = sqrt(2)*D[17]/2 + W[8]/2;
	M[76] = D[13];
	M[77] = sqrt(2)*D[15]/2 - W[6]/2;
	M[78] = sqrt(2)*D[16]/2 - W[7]/2;
	M[79] = sqrt(2)*D[15]/2 + W[6]/2;
	M[80] = D[14];

  return 0;
}

int truesdell_tangent_outer(const double * const S, double * const M)
{
	M[0] = S[0];
	M[1] = sqrt(2)*S[5];
	M[2] = sqrt(2)*S[4];
	M[3] = 0;
	M[4] = -S[0];
	M[5] = 0;
	M[6] = 0;
	M[7] = 0;
	M[8] = -S[0];
	M[9] = 0;
	M[10] = S[1];
	M[11] = sqrt(2)*S[3]/2;
	M[12] = S[0];
	M[13] = 0;
	M[14] = sqrt(2)*S[4]/2;
	M[15] = 0;
	M[16] = 0;
	M[17] = -sqrt(2)*S[5]/2;
	M[18] = 0;
	M[19] = sqrt(2)*S[3]/2;
	M[20] = S[2];
	M[21] = 0;
	M[22] = -sqrt(2)*S[4]/2;
	M[23] = 0;
	M[24] = S[0];
	M[25] = sqrt(2)*S[5]/2;
	M[26] = 0;
	M[27] = 0;
	M[28] = S[1];
	M[29] = sqrt(2)*S[3]/2;
	M[30] = S[0];
	M[31] = 0;
	M[32] = sqrt(2)*S[4]/2;
	M[33] = 0;
	M[34] = 0;
	M[35] = -sqrt(2)*S[5]/2;
	M[36] = -S[1];
	M[37] = 0;
	M[38] = 0;
	M[39] = sqrt(2)*S[5];
	M[40] = S[1];
	M[41] = sqrt(2)*S[3];
	M[42] = 0;
	M[43] = 0;
	M[44] = -S[1];
	M[45] = -sqrt(2)*S[3]/2;
	M[46] = 0;
	M[47] = 0;
	M[48] = sqrt(2)*S[4]/2;
	M[49] = 0;
	M[50] = S[2];
	M[51] = sqrt(2)*S[5]/2;
	M[52] = S[1];
	M[53] = 0;
	M[54] = 0;
	M[55] = sqrt(2)*S[3]/2;
	M[56] = S[2];
	M[57] = 0;
	M[58] = -sqrt(2)*S[4]/2;
	M[59] = 0;
	M[60] = S[0];
	M[61] = sqrt(2)*S[5]/2;
	M[62] = 0;
	M[63] = -sqrt(2)*S[3]/2;
	M[64] = 0;
	M[65] = 0;
	M[66] = sqrt(2)*S[4]/2;
	M[67] = 0;
	M[68] = S[2];
	M[69] = sqrt(2)*S[5]/2;
	M[70] = S[1];
	M[71] = 0;
	M[72] = -S[2];
	M[73] = 0;
	M[74] = 0;
	M[75] = 0;
	M[76] = -S[2];
	M[77] = 0;
	M[78] = sqrt(2)*S[4];
	M[79] = sqrt(2)*S[3];
	M[80] = S[2];
  
  return 0;
}

int full2skew(const double * const A, double * const M)
{
	M[0] = -A[5];
	M[1] = A[2];
	M[2] = -A[1];
	M[3] = -A[41];
	M[4] = A[38];
	M[5] = -A[37];
	M[6] = -A[77];
	M[7] = A[74];
	M[8] = -A[73];
	M[9] = -sqrt(2)*A[50];
	M[10] = sqrt(2)*A[47];
	M[11] = -sqrt(2)*A[46];
	M[12] = -sqrt(2)*A[23];
	M[13] = sqrt(2)*A[20];
	M[14] = -sqrt(2)*A[19];
	M[15] = -sqrt(2)*A[14];
	M[16] = sqrt(2)*A[11];
	M[17] = -sqrt(2)*A[10];

  return 0;
}

int skew2full(const double * const M, double * const A)
{
	A[0] = 0;
	A[1] = -M[2];
	A[2] = M[1];
	A[3] = M[2];
	A[4] = 0;
	A[5] = -M[0];
	A[6] = -M[1];
	A[7] = M[0];
	A[8] = 0;
	A[9] = 0;
	A[10] = -sqrt(2)*M[17]/2;
	A[11] = sqrt(2)*M[16]/2;
	A[12] = sqrt(2)*M[17]/2;
	A[13] = 0;
	A[14] = -sqrt(2)*M[15]/2;
	A[15] = -sqrt(2)*M[16]/2;
	A[16] = sqrt(2)*M[15]/2;
	A[17] = 0;
	A[18] = 0;
	A[19] = -sqrt(2)*M[14]/2;
	A[20] = sqrt(2)*M[13]/2;
	A[21] = sqrt(2)*M[14]/2;
	A[22] = 0;
	A[23] = -sqrt(2)*M[12]/2;
	A[24] = -sqrt(2)*M[13]/2;
	A[25] = sqrt(2)*M[12]/2;
	A[26] = 0;
	A[27] = 0;
	A[28] = -sqrt(2)*M[17]/2;
	A[29] = sqrt(2)*M[16]/2;
	A[30] = sqrt(2)*M[17]/2;
	A[31] = 0;
	A[32] = -sqrt(2)*M[15]/2;
	A[33] = -sqrt(2)*M[16]/2;
	A[34] = sqrt(2)*M[15]/2;
	A[35] = 0;
	A[36] = 0;
	A[37] = -M[5];
	A[38] = M[4];
	A[39] = M[5];
	A[40] = 0;
	A[41] = -M[3];
	A[42] = -M[4];
	A[43] = M[3];
	A[44] = 0;
	A[45] = 0;
	A[46] = -sqrt(2)*M[11]/2;
	A[47] = sqrt(2)*M[10]/2;
	A[48] = sqrt(2)*M[11]/2;
	A[49] = 0;
	A[50] = -sqrt(2)*M[9]/2;
	A[51] = -sqrt(2)*M[10]/2;
	A[52] = sqrt(2)*M[9]/2;
	A[53] = 0;
	A[54] = 0;
	A[55] = -sqrt(2)*M[14]/2;
	A[56] = sqrt(2)*M[13]/2;
	A[57] = sqrt(2)*M[14]/2;
	A[58] = 0;
	A[59] = -sqrt(2)*M[12]/2;
	A[60] = -sqrt(2)*M[13]/2;
	A[61] = sqrt(2)*M[12]/2;
	A[62] = 0;
	A[63] = 0;
	A[64] = -sqrt(2)*M[11]/2;
	A[65] = sqrt(2)*M[10]/2;
	A[66] = sqrt(2)*M[11]/2;
	A[67] = 0;
	A[68] = -sqrt(2)*M[9]/2;
	A[69] = -sqrt(2)*M[10]/2;
	A[70] = sqrt(2)*M[9]/2;
	A[71] = 0;
	A[72] = 0;
	A[73] = -M[8];
	A[74] = M[7];
	A[75] = M[8];
	A[76] = 0;
	A[77] = -M[6];
	A[78] = -M[7];
	A[79] = M[6];
	A[80] = 0;

  return 0;
}

int full2mandel(const double * const A, double * const M)
{
	M[0] = A[0];
	M[1] = A[4];
	M[2] = A[8];
	M[3] = sqrt(2)*A[5];
	M[4] = sqrt(2)*A[2];
	M[5] = sqrt(2)*A[1];
	M[6] = A[36];
	M[7] = A[40];
	M[8] = A[44];
	M[9] = sqrt(2)*A[41];
	M[10] = sqrt(2)*A[38];
	M[11] = sqrt(2)*A[37];
	M[12] = A[72];
	M[13] = A[76];
	M[14] = A[80];
	M[15] = sqrt(2)*A[77];
	M[16] = sqrt(2)*A[74];
	M[17] = sqrt(2)*A[73];
	M[18] = sqrt(2)*A[45];
	M[19] = sqrt(2)*A[49];
	M[20] = sqrt(2)*A[53];
	M[21] = 2*A[50];
	M[22] = 2*A[47];
	M[23] = 2*A[46];
	M[24] = sqrt(2)*A[18];
	M[25] = sqrt(2)*A[22];
	M[26] = sqrt(2)*A[26];
	M[27] = 2*A[23];
	M[28] = 2*A[20];
	M[29] = 2*A[19];
	M[30] = sqrt(2)*A[9];
	M[31] = sqrt(2)*A[13];
	M[32] = sqrt(2)*A[17];
	M[33] = 2*A[14];
	M[34] = 2*A[11];
	M[35] = 2*A[10];

  return 0;
}

int mandel2full(const double * const M, double * const A)
{
	A[0] = M[0];
	A[1] = sqrt(2)*M[5]/2;
	A[2] = sqrt(2)*M[4]/2;
	A[3] = sqrt(2)*M[5]/2;
	A[4] = M[1];
	A[5] = sqrt(2)*M[3]/2;
	A[6] = sqrt(2)*M[4]/2;
	A[7] = sqrt(2)*M[3]/2;
	A[8] = M[2];
	A[9] = sqrt(2)*M[30]/2;
	A[10] = M[35]/2;
	A[11] = M[34]/2;
	A[12] = M[35]/2;
	A[13] = sqrt(2)*M[31]/2;
	A[14] = M[33]/2;
	A[15] = M[34]/2;
	A[16] = M[33]/2;
	A[17] = sqrt(2)*M[32]/2;
	A[18] = sqrt(2)*M[24]/2;
	A[19] = M[29]/2;
	A[20] = M[28]/2;
	A[21] = M[29]/2;
	A[22] = sqrt(2)*M[25]/2;
	A[23] = M[27]/2;
	A[24] = M[28]/2;
	A[25] = M[27]/2;
	A[26] = sqrt(2)*M[26]/2;
	A[27] = sqrt(2)*M[30]/2;
	A[28] = M[35]/2;
	A[29] = M[34]/2;
	A[30] = M[35]/2;
	A[31] = sqrt(2)*M[31]/2;
	A[32] = M[33]/2;
	A[33] = M[34]/2;
	A[34] = M[33]/2;
	A[35] = sqrt(2)*M[32]/2;
	A[36] = M[6];
	A[37] = sqrt(2)*M[11]/2;
	A[38] = sqrt(2)*M[10]/2;
	A[39] = sqrt(2)*M[11]/2;
	A[40] = M[7];
	A[41] = sqrt(2)*M[9]/2;
	A[42] = sqrt(2)*M[10]/2;
	A[43] = sqrt(2)*M[9]/2;
	A[44] = M[8];
	A[45] = sqrt(2)*M[18]/2;
	A[46] = M[23]/2;
	A[47] = M[22]/2;
	A[48] = M[23]/2;
	A[49] = sqrt(2)*M[19]/2;
	A[50] = M[21]/2;
	A[51] = M[22]/2;
	A[52] = M[21]/2;
	A[53] = sqrt(2)*M[20]/2;
	A[54] = sqrt(2)*M[24]/2;
	A[55] = M[29]/2;
	A[56] = M[28]/2;
	A[57] = M[29]/2;
	A[58] = sqrt(2)*M[25]/2;
	A[59] = M[27]/2;
	A[60] = M[28]/2;
	A[61] = M[27]/2;
	A[62] = sqrt(2)*M[26]/2;
	A[63] = sqrt(2)*M[18]/2;
	A[64] = M[23]/2;
	A[65] = M[22]/2;
	A[66] = M[23]/2;
	A[67] = sqrt(2)*M[19]/2;
	A[68] = M[21]/2;
	A[69] = M[22]/2;
	A[70] = M[21]/2;
	A[71] = sqrt(2)*M[20]/2;
	A[72] = M[12];
	A[73] = sqrt(2)*M[17]/2;
	A[74] = sqrt(2)*M[16]/2;
	A[75] = sqrt(2)*M[17]/2;
	A[76] = M[13];
	A[77] = sqrt(2)*M[15]/2;
	A[78] = sqrt(2)*M[16]/2;
	A[79] = sqrt(2)*M[15]/2;
	A[80] = M[14];
  
  return 0;
}

int truesdell_update_sym(const double * const D, const double * const W,
                         const double * const Sn, const double * const So,
                         double * const Snp1)
{
  double rhs_mandel[6];
  double rhs_full[9];
  double mat[81];

  truesdell_rhs(D, W, Sn, So, rhs_mandel);
  usym(rhs_mandel, rhs_full);

  truesdell_mat(D, W, mat);

  int ier = solve_mat(mat, 9, rhs_full);
  sym(rhs_full, rhs_mandel);

  add_vec(Sn, rhs_mandel, 6, Snp1);

  return ier;
}

int truesdell_mat(const double * const D, const double * const W,
                  double * const M)
{
	M[0] = -D[0] + D[1] + D[2] + 1;
	M[1] = -sqrt(2)*D[5]/2 + W[2];
	M[2] = -sqrt(2)*D[4]/2 - W[1];
	M[3] = -sqrt(2)*D[5]/2 + W[2];
	M[4] = 0;
	M[5] = 0;
	M[6] = -sqrt(2)*D[4]/2 - W[1];
	M[7] = 0;
	M[8] = 0;
	M[9] = -sqrt(2)*D[5]/2 - W[2];
	M[10] = D[2] + 1;
	M[11] = -sqrt(2)*D[3]/2 + W[0];
	M[12] = 0;
	M[13] = -sqrt(2)*D[5]/2 + W[2];
	M[14] = 0;
	M[15] = 0;
	M[16] = -sqrt(2)*D[4]/2 - W[1];
	M[17] = 0;
	M[18] = -sqrt(2)*D[4]/2 + W[1];
	M[19] = -sqrt(2)*D[3]/2 - W[0];
	M[20] = D[1] + 1;
	M[21] = 0;
	M[22] = 0;
	M[23] = -sqrt(2)*D[5]/2 + W[2];
	M[24] = 0;
	M[25] = 0;
	M[26] = -sqrt(2)*D[4]/2 - W[1];
	M[27] = -sqrt(2)*D[5]/2 - W[2];
	M[28] = 0;
	M[29] = 0;
	M[30] = D[2] + 1;
	M[31] = -sqrt(2)*D[5]/2 + W[2];
	M[32] = -sqrt(2)*D[4]/2 - W[1];
	M[33] = -sqrt(2)*D[3]/2 + W[0];
	M[34] = 0;
	M[35] = 0;
	M[36] = 0;
	M[37] = -sqrt(2)*D[5]/2 - W[2];
	M[38] = 0;
	M[39] = -sqrt(2)*D[5]/2 - W[2];
	M[40] = D[0] - D[1] + D[2] + 1;
	M[41] = -sqrt(2)*D[3]/2 + W[0];
	M[42] = 0;
	M[43] = -sqrt(2)*D[3]/2 + W[0];
	M[44] = 0;
	M[45] = 0;
	M[46] = 0;
	M[47] = -sqrt(2)*D[5]/2 - W[2];
	M[48] = -sqrt(2)*D[4]/2 + W[1];
	M[49] = -sqrt(2)*D[3]/2 - W[0];
	M[50] = D[0] + 1;
	M[51] = 0;
	M[52] = 0;
	M[53] = -sqrt(2)*D[3]/2 + W[0];
	M[54] = -sqrt(2)*D[4]/2 + W[1];
	M[55] = 0;
	M[56] = 0;
	M[57] = -sqrt(2)*D[3]/2 - W[0];
	M[58] = 0;
	M[59] = 0;
	M[60] = D[1] + 1;
	M[61] = -sqrt(2)*D[5]/2 + W[2];
	M[62] = -sqrt(2)*D[4]/2 - W[1];
	M[63] = 0;
	M[64] = -sqrt(2)*D[4]/2 + W[1];
	M[65] = 0;
	M[66] = 0;
	M[67] = -sqrt(2)*D[3]/2 - W[0];
	M[68] = 0;
	M[69] = -sqrt(2)*D[5]/2 - W[2];
	M[70] = D[0] + 1;
	M[71] = -sqrt(2)*D[3]/2 + W[0];
	M[72] = 0;
	M[73] = 0;
	M[74] = -sqrt(2)*D[4]/2 + W[1];
	M[75] = 0;
	M[76] = 0;
	M[77] = -sqrt(2)*D[3]/2 - W[0];
	M[78] = -sqrt(2)*D[4]/2 + W[1];
	M[79] = -sqrt(2)*D[3]/2 - W[0];
	M[80] = D[0] + D[1] - D[2] + 1;

  return 0;
}

int truesdell_rhs(const double * const D, const double * const W,
                  const double * const Sn, const double * const So,
                  double * const St)
{
	St[0] = D[0]*Sn[0] - D[1]*Sn[0] - D[2]*Sn[0] + D[4]*Sn[4] + D[5]*Sn[5] + sqrt(2)*Sn[4]*W[1] - sqrt(2)*Sn[5]*W[2] + So[0];
	St[1] = -D[0]*Sn[1] + D[1]*Sn[1] - D[2]*Sn[1] + D[3]*Sn[3] + D[5]*Sn[5] - sqrt(2)*Sn[3]*W[0] + sqrt(2)*Sn[5]*W[2] + So[1];
	St[2] = -D[0]*Sn[2] - D[1]*Sn[2] + D[2]*Sn[2] + D[3]*Sn[3] + D[4]*Sn[4] + sqrt(2)*Sn[3]*W[0] - sqrt(2)*Sn[4]*W[1] + So[2];
	St[3] = sqrt(2)*(-sqrt(2)*D[0]*Sn[3]/2 + sqrt(2)*D[3]*Sn[1]/2 + sqrt(2)*D[3]*Sn[2]/2 + D[4]*Sn[5]/2 + D[5]*Sn[4]/2 + Sn[1]*W[0] - Sn[2]*W[0] + sqrt(2)*Sn[4]*W[2]/2 - sqrt(2)*Sn[5]*W[1]/2 + sqrt(2)*So[3]/2);
	St[4] = sqrt(2)*(-sqrt(2)*D[1]*Sn[4]/2 + D[3]*Sn[5]/2 + sqrt(2)*D[4]*Sn[0]/2 + sqrt(2)*D[4]*Sn[2]/2 + D[5]*Sn[3]/2 - Sn[0]*W[1] + Sn[2]*W[1] - sqrt(2)*Sn[3]*W[2]/2 + sqrt(2)*Sn[5]*W[0]/2 + sqrt(2)*So[4]/2);
	St[5] = sqrt(2)*(-sqrt(2)*D[2]*Sn[5]/2 + D[3]*Sn[4]/2 + D[4]*Sn[3]/2 + sqrt(2)*D[5]*Sn[0]/2 + sqrt(2)*D[5]*Sn[1]/2 + Sn[0]*W[2] - Sn[1]*W[2] + sqrt(2)*Sn[3]*W[1]/2 - sqrt(2)*Sn[4]*W[0]/2 + sqrt(2)*So[5]/2);

  return 0;
}

int sym(const double * const A, double * const v)
{
  v[0] = A[0];
  v[1] = A[4];
  v[2] = A[8];
  v[3] = sqrt(2.0) * A[5];
  v[4] = sqrt(2.0) * A[2];
  v[5] = sqrt(2.0) * A[1];

  return 0;
}

int usym(const double * const v, double * const A)
{
  A[0] = v[0];
  A[1] = v[5] / sqrt(2.0);
  A[2] = v[4] / sqrt(2.0);
  A[3] = v[5] / sqrt(2.0);
  A[4] = v[1];
  A[5] = v[3] / sqrt(2.0);
  A[6] = v[4] / sqrt(2.0);
  A[7] = v[3] / sqrt(2.0);
  A[8] = v[2];

  return 0;
}

int minus_vec(double * const a, int n)
{
  for (int i=0; i<n; i++) {
    a[i] = -a[i];
  }

  return 0;
}

int add_vec(const double * const a, const double * const b, int n, double * const c)
{
  for (int i=0; i<n; i++) {
    c[i] = a[i] + b[i];
  }

  return 0;
}

int sub_vec(const double * const a, const double * const b, int n, double * const c)
{
  for (int i=0; i<n; i++) {
    c[i] = a[i] - b[i];
  }
  
  return 0;
}

double dot_vec(const double * const a, const double * const b, int n)
{
  double sum = 0.0;
  for (int i=0; i<n; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

double norm2_vec(const double * const a, int n)
{
  return sqrt(dot_vec(a, a, n));
}

int normalize_vec(double * const a, int n)
{
  double nv = norm2_vec(a, n);
  if (fabs(nv) < std::numeric_limits<double>::epsilon()) {
    std::fill(a, a+n, 0.0);
    return 0;
  }
  for (int i=0; i<n; i++) {
    a[i] /= nv;
  }
  return 0;
}

int dev_vec(double * const a)
{
  double tr = (a[0] + a[1] + a[2]) / 3.0;
  for (int i=0; i<3; i++) {
    a[i] -= tr;
  }

  return 0;
}

int outer_vec(const double * const a, int na, const double * const b, int nb, double * const C)
{
  for (int i=0; i < na; i++) {
    for (int j=0; j < nb; j++) {
      C[CINDEX(i,j,nb)] = a[i] * b[j];
    }
  }

  return 0;
}

// Problem?
int outer_update(const double * const a, int na, const double * const b, 
                 int nb, double * const C)
{
  dger_(nb, na, 1.0, b, 1, a, 1, C, nb);

  return 0;
}

int outer_update_minus(const double * const a, int na, const double * const b, 
                       int nb, double * const C)
{
  dger_(nb, na, -1.0, b, 1, a, 1, C, nb);

  return 0;
}

int mat_vec(const double * const A, int m, const double * const b, int n, 
            double * const c)
{
  dgemv_("T", n, m, 1.0, A, n, b, 1, 0.0, c, 1);

  return 0;
}

int mat_vec_trans(const double * const A, int m, const double * const b, int n, 
            double * const c)
{
  dgemv_("N", m, n, 1.0, A, m, b, 1, 0.0, c, 1);

  return 0;
}

int invert_mat(double * const A, int n)
{
  int * ipiv = new int[n + 1];
  int lwork = n * n;
  double * work = new double[lwork];
  int info;

  dgetrf_(n, n, A, n, ipiv, info);
  if (info > 0) {
    delete [] ipiv;
    delete [] work;
    return LINALG_FAILURE;
  }

  dgetri_(n, A, n, ipiv, work, lwork, info);

  delete [] ipiv;
  delete [] work;

  if (info > 0) return LINALG_FAILURE;

  return 0;
}

int mat_mat(int m, int n, int k, const double * const A,
            const double * const B, double * const C)
{
  dgemm_("N", "N", n, m, k, 1.0, B, n, A, k, 0.0, C, n);

  return 0;
}

int solve_mat(const double * const A, int n, double * const x)
{
  int info;
  int * ipiv = new int [n];
  double * B = new double [n*n];
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      B[CINDEX(i,j,n)] = A[CINDEX(j,i,n)];
    }
  }
  
  dgesv_(n, 1, B, n, ipiv, x, n, info);

  delete [] ipiv;
  delete [] B;

  if (info > 0) return LINALG_FAILURE;
  
  return 0;
}

/*
 *  No error checking in this function, as it is assumed to be non-critical
 */
double condition(const double * const A, int n)
{
  // Setup
  int info;
  int * ipiv = new int [n];
  double * x = new double [n];
  double * B = new double [n*n];
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      B[CINDEX(i,j,n)] = A[CINDEX(j,i,n)];
    }
  }
  
  // 1 norm
  double anorm = 0.0;
  double csum;
  for (int i=0; i<n; i++) {
    csum = 0.0;
    for (int j=0; j<n; j++) {
      csum += fabs(A[CINDEX(i,j,n)]);
    }
    if (csum > anorm) anorm = csum;
  }

  // Solve
  std::fill(x, x+n, 0.0);
  dgesv_(n, 1, B, n, ipiv, x, n, info);
  delete [] ipiv;

  double * work = new double [4*n];
  int * iwork = new int [n];
  double rcond;
  dgecon_("1", n, B, n, anorm, rcond, work, iwork, info);

  delete [] B;
  delete [] work;
  delete [] iwork;

  return 1.0 / rcond;
}

double polyval(const double * const poly, const int n, double x)
{
  double res = poly[0];
  for (int i=1; i < n; i++) {
    res = res * x + poly[i];
  }
  return res;
}

int gcd(int a, int b)
{
  if (a == 0) return b;
  return gcd(b % a, a);
}

int common_gcd(std::vector<int> in)
{
  int c = in[0];
  for (int i = 1; i < in.size(); i++) {
    c = gcd(c, in[i]);
  }
  return c;
}

std::vector<int> reduce_gcd(std::vector<int> in)
{
  std::vector<int> out(in);
  int f = common_gcd(out);
  for (int i=0; i<out.size(); i++) {
    out[i] /= abs(f);
  }

  return out;
}

double convert_angle(double a, std::string angles)
{
  if (std::string("radians").compare(angles) == 0) {
    return a;
  }
  else if (std::string("degrees").compare(angles) == 0) {
    return a / 360.0 * 2.0 * M_PI;
  }
  else {
    throw std::invalid_argument("Angle type must be radians or degrees.");
  }
}

double cast_angle(double a, std::string angles)
{
  if (std::string("radians").compare(angles) == 0) {
    return a;
  }
  else if (std::string("degrees").compare(angles) == 0) {
    return a / (2.0 * M_PI) * 360.0;
  }
  else {
    throw std::invalid_argument("Angle type must be radians or degrees.");
  }
}

} // namespace neml
