#include "math/nemlmath.h"

#include "nemlerror.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

namespace neml {

void SymSymR4SkewmSkewSymR4SymR4(const double * const M, const double * const W, double * const SS)
{
	SS[0] = sqrt(2)*(-M[24]*W[1] + M[30]*W[2]);
	SS[1] = sqrt(2)*(-M[25]*W[1] + M[31]*W[2]);
	SS[2] = sqrt(2)*(-M[26]*W[1] + M[32]*W[2]);
	SS[3] = sqrt(2)*(-M[27]*W[1] + M[33]*W[2]);
	SS[4] = sqrt(2)*(-M[28]*W[1] + M[34]*W[2]);
	SS[5] = sqrt(2)*(-M[29]*W[1] + M[35]*W[2]);
	SS[6] = sqrt(2)*(M[18]*W[0] - M[30]*W[2]);
	SS[7] = sqrt(2)*(M[19]*W[0] - M[31]*W[2]);
	SS[8] = sqrt(2)*(M[20]*W[0] - M[32]*W[2]);
	SS[9] = sqrt(2)*(M[21]*W[0] - M[33]*W[2]);
	SS[10] = sqrt(2)*(M[22]*W[0] - M[34]*W[2]);
	SS[11] = sqrt(2)*(M[23]*W[0] - M[35]*W[2]);
	SS[12] = sqrt(2)*(-M[18]*W[0] + M[24]*W[1]);
	SS[13] = sqrt(2)*(-M[19]*W[0] + M[25]*W[1]);
	SS[14] = sqrt(2)*(-M[20]*W[0] + M[26]*W[1]);
	SS[15] = sqrt(2)*(-M[21]*W[0] + M[27]*W[1]);
	SS[16] = sqrt(2)*(-M[22]*W[0] + M[28]*W[1]);
	SS[17] = sqrt(2)*(-M[23]*W[0] + M[29]*W[1]);
	SS[18] = sqrt(2)*M[12]*W[0] - M[24]*W[2] + M[30]*W[1] - sqrt(2)*M[6]*W[0];
	SS[19] = sqrt(2)*M[13]*W[0] - M[25]*W[2] + M[31]*W[1] - sqrt(2)*M[7]*W[0];
	SS[20] = sqrt(2)*M[14]*W[0] - M[26]*W[2] + M[32]*W[1] - sqrt(2)*M[8]*W[0];
	SS[21] = sqrt(2)*M[15]*W[0] - M[27]*W[2] + M[33]*W[1] - sqrt(2)*M[9]*W[0];
	SS[22] = -sqrt(2)*M[10]*W[0] + sqrt(2)*M[16]*W[0] - M[28]*W[2] + M[34]*W[1];
	SS[23] = -sqrt(2)*M[11]*W[0] + sqrt(2)*M[17]*W[0] - M[29]*W[2] + M[35]*W[1];
	SS[24] = sqrt(2)*M[0]*W[1] - sqrt(2)*M[12]*W[1] + M[18]*W[2] - M[30]*W[0];
	SS[25] = -sqrt(2)*M[13]*W[1] + M[19]*W[2] + sqrt(2)*M[1]*W[1] - M[31]*W[0];
	SS[26] = -sqrt(2)*M[14]*W[1] + M[20]*W[2] + sqrt(2)*M[2]*W[1] - M[32]*W[0];
	SS[27] = -sqrt(2)*M[15]*W[1] + M[21]*W[2] - M[33]*W[0] + sqrt(2)*M[3]*W[1];
	SS[28] = -sqrt(2)*M[16]*W[1] + M[22]*W[2] - M[34]*W[0] + sqrt(2)*M[4]*W[1];
	SS[29] = -sqrt(2)*M[17]*W[1] + M[23]*W[2] - M[35]*W[0] + sqrt(2)*M[5]*W[1];
	SS[30] = -sqrt(2)*M[0]*W[2] - M[18]*W[1] + M[24]*W[0] + sqrt(2)*M[6]*W[2];
	SS[31] = -M[19]*W[1] - sqrt(2)*M[1]*W[2] + M[25]*W[0] + sqrt(2)*M[7]*W[2];
	SS[32] = -M[20]*W[1] + M[26]*W[0] - sqrt(2)*M[2]*W[2] + sqrt(2)*M[8]*W[2];
	SS[33] = -M[21]*W[1] + M[27]*W[0] - sqrt(2)*M[3]*W[2] + sqrt(2)*M[9]*W[2];
	SS[34] = sqrt(2)*M[10]*W[2] - M[22]*W[1] + M[28]*W[0] - sqrt(2)*M[4]*W[2];
	SS[35] = sqrt(2)*M[11]*W[2] - M[23]*W[1] + M[29]*W[0] - sqrt(2)*M[5]*W[2];
}

void SymSkewR4SymmSkewSymR4SymR4(const double * const D, const double * const M, double * const SS)
{
	SS[0] = sqrt(2)*(-D[4]*M[6] + D[5]*M[12]);
	SS[1] = sqrt(2)*(-D[4]*M[7] + D[5]*M[13]);
	SS[2] = sqrt(2)*(-D[4]*M[8] + D[5]*M[14]);
	SS[3] = sqrt(2)*(-D[4]*M[9] + D[5]*M[15]);
	SS[4] = sqrt(2)*(-D[4]*M[10] + D[5]*M[16]);
	SS[5] = sqrt(2)*(-D[4]*M[11] + D[5]*M[17]);
	SS[6] = sqrt(2)*(D[3]*M[0] - D[5]*M[12]);
	SS[7] = sqrt(2)*(D[3]*M[1] - D[5]*M[13]);
	SS[8] = sqrt(2)*(D[3]*M[2] - D[5]*M[14]);
	SS[9] = sqrt(2)*(D[3]*M[3] - D[5]*M[15]);
	SS[10] = sqrt(2)*(D[3]*M[4] - D[5]*M[16]);
	SS[11] = sqrt(2)*(D[3]*M[5] - D[5]*M[17]);
	SS[12] = sqrt(2)*(-D[3]*M[0] + D[4]*M[6]);
	SS[13] = sqrt(2)*(-D[3]*M[1] + D[4]*M[7]);
	SS[14] = sqrt(2)*(-D[3]*M[2] + D[4]*M[8]);
	SS[15] = sqrt(2)*(-D[3]*M[3] + D[4]*M[9]);
	SS[16] = sqrt(2)*(-D[3]*M[4] + D[4]*M[10]);
	SS[17] = sqrt(2)*(-D[3]*M[5] + D[4]*M[11]);
	SS[18] = -sqrt(2)*D[1]*M[0] + sqrt(2)*D[2]*M[0] - D[4]*M[12] + D[5]*M[6];
	SS[19] = -sqrt(2)*D[1]*M[1] + sqrt(2)*D[2]*M[1] - D[4]*M[13] + D[5]*M[7];
	SS[20] = -sqrt(2)*D[1]*M[2] + sqrt(2)*D[2]*M[2] - D[4]*M[14] + D[5]*M[8];
	SS[21] = -sqrt(2)*D[1]*M[3] + sqrt(2)*D[2]*M[3] - D[4]*M[15] + D[5]*M[9];
	SS[22] = -sqrt(2)*D[1]*M[4] + sqrt(2)*D[2]*M[4] - D[4]*M[16] + D[5]*M[10];
	SS[23] = -sqrt(2)*D[1]*M[5] + sqrt(2)*D[2]*M[5] - D[4]*M[17] + D[5]*M[11];
	SS[24] = sqrt(2)*D[0]*M[6] - sqrt(2)*D[2]*M[6] + D[3]*M[12] - D[5]*M[0];
	SS[25] = sqrt(2)*D[0]*M[7] - sqrt(2)*D[2]*M[7] + D[3]*M[13] - D[5]*M[1];
	SS[26] = sqrt(2)*D[0]*M[8] - sqrt(2)*D[2]*M[8] + D[3]*M[14] - D[5]*M[2];
	SS[27] = sqrt(2)*D[0]*M[9] - sqrt(2)*D[2]*M[9] + D[3]*M[15] - D[5]*M[3];
	SS[28] = sqrt(2)*D[0]*M[10] - sqrt(2)*D[2]*M[10] + D[3]*M[16] - D[5]*M[4];
	SS[29] = sqrt(2)*D[0]*M[11] - sqrt(2)*D[2]*M[11] + D[3]*M[17] - D[5]*M[5];
	SS[30] = -sqrt(2)*D[0]*M[12] + sqrt(2)*D[1]*M[12] - D[3]*M[6] + D[4]*M[0];
	SS[31] = -sqrt(2)*D[0]*M[13] + sqrt(2)*D[1]*M[13] - D[3]*M[7] + D[4]*M[1];
	SS[32] = -sqrt(2)*D[0]*M[14] + sqrt(2)*D[1]*M[14] - D[3]*M[8] + D[4]*M[2];
	SS[33] = -sqrt(2)*D[0]*M[15] + sqrt(2)*D[1]*M[15] - D[3]*M[9] + D[4]*M[3];
	SS[34] = -sqrt(2)*D[0]*M[16] + sqrt(2)*D[1]*M[16] - D[3]*M[10] + D[4]*M[4];
	SS[35] = -sqrt(2)*D[0]*M[17] + sqrt(2)*D[1]*M[17] - D[3]*M[11] + D[4]*M[5];
}

void SpecialSymSymR4Sym(const double * const D, const double * const M, double * const SW)
{
	SW[0] = -sqrt(2)*D[1]*M[3]/2 + sqrt(2)*D[2]*M[3]/2 + sqrt(2)*D[3]*M[1]/2 - sqrt(2)*D[3]*M[2]/2 + D[4]*M[5]/2 - D[5]*M[4]/2;
	SW[1] = sqrt(2)*D[0]*M[4]/2 - sqrt(2)*D[2]*M[4]/2 - D[3]*M[5]/2 - sqrt(2)*D[4]*M[0]/2 + sqrt(2)*D[4]*M[2]/2 + D[5]*M[3]/2;
	SW[2] = -sqrt(2)*D[0]*M[5]/2 + sqrt(2)*D[1]*M[5]/2 + D[3]*M[4]/2 - D[4]*M[3]/2 + sqrt(2)*D[5]*M[0]/2 - sqrt(2)*D[5]*M[1]/2;
	SW[3] = -sqrt(2)*D[1]*M[9]/2 + sqrt(2)*D[2]*M[9]/2 + sqrt(2)*D[3]*M[7]/2 - sqrt(2)*D[3]*M[8]/2 + D[4]*M[11]/2 - D[5]*M[10]/2;
	SW[4] = sqrt(2)*D[0]*M[10]/2 - sqrt(2)*D[2]*M[10]/2 - D[3]*M[11]/2 - sqrt(2)*D[4]*M[6]/2 + sqrt(2)*D[4]*M[8]/2 + D[5]*M[9]/2;
	SW[5] = -sqrt(2)*D[0]*M[11]/2 + sqrt(2)*D[1]*M[11]/2 + D[3]*M[10]/2 - D[4]*M[9]/2 + sqrt(2)*D[5]*M[6]/2 - sqrt(2)*D[5]*M[7]/2;
	SW[6] = -sqrt(2)*D[1]*M[15]/2 + sqrt(2)*D[2]*M[15]/2 + sqrt(2)*D[3]*M[13]/2 - sqrt(2)*D[3]*M[14]/2 + D[4]*M[17]/2 - D[5]*M[16]/2;
	SW[7] = sqrt(2)*D[0]*M[16]/2 - sqrt(2)*D[2]*M[16]/2 - D[3]*M[17]/2 - sqrt(2)*D[4]*M[12]/2 + sqrt(2)*D[4]*M[14]/2 + D[5]*M[15]/2;
	SW[8] = -sqrt(2)*D[0]*M[17]/2 + sqrt(2)*D[1]*M[17]/2 + D[3]*M[16]/2 - D[4]*M[15]/2 + sqrt(2)*D[5]*M[12]/2 - sqrt(2)*D[5]*M[13]/2;
	SW[9] = -sqrt(2)*D[1]*M[21]/2 + sqrt(2)*D[2]*M[21]/2 + sqrt(2)*D[3]*M[19]/2 - sqrt(2)*D[3]*M[20]/2 + D[4]*M[23]/2 - D[5]*M[22]/2;
	SW[10] = sqrt(2)*D[0]*M[22]/2 - sqrt(2)*D[2]*M[22]/2 - D[3]*M[23]/2 - sqrt(2)*D[4]*M[18]/2 + sqrt(2)*D[4]*M[20]/2 + D[5]*M[21]/2;
	SW[11] = -sqrt(2)*D[0]*M[23]/2 + sqrt(2)*D[1]*M[23]/2 + D[3]*M[22]/2 - D[4]*M[21]/2 + sqrt(2)*D[5]*M[18]/2 - sqrt(2)*D[5]*M[19]/2;
	SW[12] = -sqrt(2)*D[1]*M[27]/2 + sqrt(2)*D[2]*M[27]/2 + sqrt(2)*D[3]*M[25]/2 - sqrt(2)*D[3]*M[26]/2 + D[4]*M[29]/2 - D[5]*M[28]/2;
	SW[13] = sqrt(2)*D[0]*M[28]/2 - sqrt(2)*D[2]*M[28]/2 - D[3]*M[29]/2 - sqrt(2)*D[4]*M[24]/2 + sqrt(2)*D[4]*M[26]/2 + D[5]*M[27]/2;
	SW[14] = -sqrt(2)*D[0]*M[29]/2 + sqrt(2)*D[1]*M[29]/2 + D[3]*M[28]/2 - D[4]*M[27]/2 + sqrt(2)*D[5]*M[24]/2 - sqrt(2)*D[5]*M[25]/2;
	SW[15] = -sqrt(2)*D[1]*M[33]/2 + sqrt(2)*D[2]*M[33]/2 + sqrt(2)*D[3]*M[31]/2 - sqrt(2)*D[3]*M[32]/2 + D[4]*M[35]/2 - D[5]*M[34]/2;
	SW[16] = sqrt(2)*D[0]*M[34]/2 - sqrt(2)*D[2]*M[34]/2 - D[3]*M[35]/2 - sqrt(2)*D[4]*M[30]/2 + sqrt(2)*D[4]*M[32]/2 + D[5]*M[33]/2;
	SW[17] = -sqrt(2)*D[0]*M[35]/2 + sqrt(2)*D[1]*M[35]/2 + D[3]*M[34]/2 - D[4]*M[33]/2 + sqrt(2)*D[5]*M[30]/2 - sqrt(2)*D[5]*M[31]/2;
}

void transform_fourth(const double * const D, const double * const W, double * const M)
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
}

void truesdell_tangent_outer(const double * const S, double * const M)
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
}

void full2skew(const double * const A, double * const M)
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
}

void skew2full(const double * const M, double * const A)
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
}

void full2wws(const double * const A, double * const M)
{
	M[0] = -A[45];
	M[1] = -A[49];
	M[2] = -A[53];
	M[3] = -sqrt(2)*A[50];
	M[4] = -sqrt(2)*A[47];
	M[5] = -sqrt(2)*A[46];
	M[6] = A[18];
	M[7] = A[22];
	M[8] = A[26];
	M[9] = sqrt(2)*A[23];
	M[10] = sqrt(2)*A[20];
	M[11] = sqrt(2)*A[19];
	M[12] = -A[9];
	M[13] = -A[13];
	M[14] = -A[17];
	M[15] = -sqrt(2)*A[14];
	M[16] = -sqrt(2)*A[11];
	M[17] = -sqrt(2)*A[10];
}

void wws2full(const double * const M, double * const A)
{
	A[0] = 0;
	A[1] = 0;
	A[2] = 0;
	A[3] = 0;
	A[4] = 0;
	A[5] = 0;
	A[6] = 0;
	A[7] = 0;
	A[8] = 0;
	A[9] = -M[12];
	A[10] = -sqrt(2)*M[17]/2;
	A[11] = -sqrt(2)*M[16]/2;
	A[12] = -sqrt(2)*M[17]/2;
	A[13] = -M[13];
	A[14] = -sqrt(2)*M[15]/2;
	A[15] = -sqrt(2)*M[16]/2;
	A[16] = -sqrt(2)*M[15]/2;
	A[17] = -M[14];
	A[18] = M[6];
	A[19] = sqrt(2)*M[11]/2;
	A[20] = sqrt(2)*M[10]/2;
	A[21] = sqrt(2)*M[11]/2;
	A[22] = M[7];
	A[23] = sqrt(2)*M[9]/2;
	A[24] = sqrt(2)*M[10]/2;
	A[25] = sqrt(2)*M[9]/2;
	A[26] = M[8];
	A[27] = M[12];
	A[28] = sqrt(2)*M[17]/2;
	A[29] = sqrt(2)*M[16]/2;
	A[30] = sqrt(2)*M[17]/2;
	A[31] = M[13];
	A[32] = sqrt(2)*M[15]/2;
	A[33] = sqrt(2)*M[16]/2;
	A[34] = sqrt(2)*M[15]/2;
	A[35] = M[14];
	A[36] = 0;
	A[37] = 0;
	A[38] = 0;
	A[39] = 0;
	A[40] = 0;
	A[41] = 0;
	A[42] = 0;
	A[43] = 0;
	A[44] = 0;
	A[45] = -M[0];
	A[46] = -sqrt(2)*M[5]/2;
	A[47] = -sqrt(2)*M[4]/2;
	A[48] = -sqrt(2)*M[5]/2;
	A[49] = -M[1];
	A[50] = -sqrt(2)*M[3]/2;
	A[51] = -sqrt(2)*M[4]/2;
	A[52] = -sqrt(2)*M[3]/2;
	A[53] = -M[2];
	A[54] = -M[6];
	A[55] = -sqrt(2)*M[11]/2;
	A[56] = -sqrt(2)*M[10]/2;
	A[57] = -sqrt(2)*M[11]/2;
	A[58] = -M[7];
	A[59] = -sqrt(2)*M[9]/2;
	A[60] = -sqrt(2)*M[10]/2;
	A[61] = -sqrt(2)*M[9]/2;
	A[62] = -M[8];
	A[63] = M[0];
	A[64] = sqrt(2)*M[5]/2;
	A[65] = sqrt(2)*M[4]/2;
	A[66] = sqrt(2)*M[5]/2;
	A[67] = M[1];
	A[68] = sqrt(2)*M[3]/2;
	A[69] = sqrt(2)*M[4]/2;
	A[70] = sqrt(2)*M[3]/2;
	A[71] = M[2];
	A[72] = 0;
	A[73] = 0;
	A[74] = 0;
	A[75] = 0;
	A[76] = 0;
	A[77] = 0;
	A[78] = 0;
	A[79] = 0;
	A[80] = 0;
}

void full2mandel(const double * const A, double * const M)
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
}

void mandel2full(const double * const M, double * const A)
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
}

void truesdell_update_sym(const double * const D, const double * const W,
                         const double * const Sn, const double * const So,
                         double * const Snp1)
{
  double rhs_mandel[6];
  double rhs_full[9];
  double mat[81];

  truesdell_rhs(D, W, Sn, So, rhs_mandel);
  usym(rhs_mandel, rhs_full);

  truesdell_mat(D, W, mat);

  solve_mat(mat, 9, rhs_full);
  sym(rhs_full, rhs_mandel);

  add_vec(Sn, rhs_mandel, 6, Snp1);
}

void truesdell_mat(const double * const D, const double * const W,
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
}

void truesdell_rhs(const double * const D, const double * const W,
                  const double * const Sn, const double * const So,
                  double * const St)
{
	St[0] = D[0]*Sn[0] - D[1]*Sn[0] - D[2]*Sn[0] + D[4]*Sn[4] + D[5]*Sn[5] + sqrt(2)*Sn[4]*W[1] - sqrt(2)*Sn[5]*W[2] + So[0];
	St[1] = -D[0]*Sn[1] + D[1]*Sn[1] - D[2]*Sn[1] + D[3]*Sn[3] + D[5]*Sn[5] - sqrt(2)*Sn[3]*W[0] + sqrt(2)*Sn[5]*W[2] + So[1];
	St[2] = -D[0]*Sn[2] - D[1]*Sn[2] + D[2]*Sn[2] + D[3]*Sn[3] + D[4]*Sn[4] + sqrt(2)*Sn[3]*W[0] - sqrt(2)*Sn[4]*W[1] + So[2];
	St[3] = sqrt(2)*(-sqrt(2)*D[0]*Sn[3]/2 + sqrt(2)*D[3]*Sn[1]/2 + sqrt(2)*D[3]*Sn[2]/2 + D[4]*Sn[5]/2 + D[5]*Sn[4]/2 + Sn[1]*W[0] - Sn[2]*W[0] + sqrt(2)*Sn[4]*W[2]/2 - sqrt(2)*Sn[5]*W[1]/2 + sqrt(2)*So[3]/2);
	St[4] = sqrt(2)*(-sqrt(2)*D[1]*Sn[4]/2 + D[3]*Sn[5]/2 + sqrt(2)*D[4]*Sn[0]/2 + sqrt(2)*D[4]*Sn[2]/2 + D[5]*Sn[3]/2 - Sn[0]*W[1] + Sn[2]*W[1] - sqrt(2)*Sn[3]*W[2]/2 + sqrt(2)*Sn[5]*W[0]/2 + sqrt(2)*So[4]/2);
	St[5] = sqrt(2)*(-sqrt(2)*D[2]*Sn[5]/2 + D[3]*Sn[4]/2 + D[4]*Sn[3]/2 + sqrt(2)*D[5]*Sn[0]/2 + sqrt(2)*D[5]*Sn[1]/2 + Sn[0]*W[2] - Sn[1]*W[2] + sqrt(2)*Sn[3]*W[1]/2 - sqrt(2)*Sn[4]*W[0]/2 + sqrt(2)*So[5]/2);
}

void sym(const double * const A, double * const v)
{
  v[0] = A[0];
  v[1] = A[4];
  v[2] = A[8];
  v[3] = sqrt(2.0) * A[5];
  v[4] = sqrt(2.0) * A[2];
  v[5] = sqrt(2.0) * A[1];
}

void usym(const double * const v, double * const A)
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
}

void skew(const double * const A, double * const v)
{
  v[0] = -A[5];
  v[1] = A[2];
  v[2] = -A[1];
}

void uskew(const double * const v, double * const A)
{
  A[0] = 0;
  A[1] = -v[2];
  A[2] = v[1];
  A[3] = v[2];
  A[4] = 0;
  A[5] = -v[0];
  A[6] = -v[1];
  A[7] = v[0];
  A[8] = 0;
}

void minus_vec(double * const a, int n)
{
  for (int i=0; i<n; i++) {
    a[i] = -a[i];
  }
}

void add_vec(const double * const a, const double * const b, int n, double * const c)
{
  for (int i=0; i<n; i++) {
    c[i] = a[i] + b[i];
  }
}

void sub_vec(const double * const a, const double * const b, int n, double * const c)
{
  for (int i=0; i<n; i++) {
    c[i] = a[i] - b[i];
  }
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

void normalize_vec(double * const a, int n)
{
  double nv = norm2_vec(a, n);
  if (fabs(nv) < std::numeric_limits<double>::epsilon()) {
    std::fill(a, a+n, 0.0);
    return;
  }
  for (int i=0; i<n; i++) {
    a[i] /= nv;
  }
}

void dev_vec(double * const a)
{
  double tr = (a[0] + a[1] + a[2]) / 3.0;
  for (int i=0; i<3; i++) {
    a[i] -= tr;
  }
}

void outer_vec(const double * const a, int na, const double * const b, int nb, double * const C)
{
  for (int i=0; i < na; i++) {
    for (int j=0; j < nb; j++) {
      C[CINDEX(i,j,nb)] = a[i] * b[j];
    }
  }
}

// Problem?
void outer_update(const double * const a, int na, const double * const b, 
                 int nb, double * const C)
{
  dger_(nb, na, 1.0, b, 1, a, 1, C, nb);
}

void outer_update_minus(const double * const a, int na, const double * const b, 
                       int nb, double * const C)
{
  dger_(nb, na, -1.0, b, 1, a, 1, C, nb);
}

void mat_vec(const double * const A, int m, const double * const b, int n, 
            double * const c)
{
  dgemv_("T", n, m, 1.0, A, n, b, 1, 0.0, c, 1);
}

void mat_vec_trans(const double * const A, int m, const double * const b, int n, 
            double * const c)
{
  dgemv_("N", m, n, 1.0, A, m, b, 1, 0.0, c, 1);
}

void invert_mat(double * const A, int n)
{
  int * ipiv = new int[n + 1];
  int lwork = n * n;
  double * work = new double[lwork];
  int info;

  dgetrf_(n, n, A, n, ipiv, info);
  if (info > 0) {
    delete [] ipiv;
    delete [] work;
    throw LinalgError("Matrix could not be inverted!");
  }

  dgetri_(n, A, n, ipiv, work, lwork, info);

  delete [] ipiv;
  delete [] work;

  if (info > 0) throw LinalgError("Matrix could not be inverted!");
}

void mat_mat(int m, int n, int k, const double * const A,
            const double * const B, double * const C)
{
  dgemm_("N", "N", n, m, k, 1.0, B, n, A, k, 0.0, C, n);
}

void mat_mat_ABT(int m, int n, int k, const double * const A,
            const double * const B, double * const C)
{
  // Provide as A_mk B_nk
  dgemm_("T", "N", n, m, k, 1.0, B, k, A, k, 0.0, C, n);
}

void solve_mat(const double * const A, int n, double * const x)
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

  if (info > 0) throw LinalgError("Matrix could not be inverted!");
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

double polyval(const std::vector<double> & poly, double x)
{
  int n = poly.size();
  double res = poly[0];
  for (int i=1; i < n; i++) {
    res = res * x + poly[i];
  }
  return res;
}

std::vector<double> poly_from_roots(const std::vector<double> & roots)
{
  std::vector<double> poly(1);
  poly[0] = 1.0;

  for (auto r = roots.begin(); r != roots.end(); ++r) {
    auto pprime = std::vector<double>(poly);
    poly.resize(poly.size()+1);
    for (int i=poly.size()-1; i > 0; i--) {
      poly[i] = poly[i-1];
    }
    poly[0] = 0.0;

    for (size_t i=0; i < pprime.size(); i++) {
      pprime[i] *= *r;
    }

    for (size_t i=0; i<(poly.size()-1); i++) {
      poly[i] -= pprime[i];
    }
  }

  std::reverse(poly.begin(), poly.end());

  return poly;
}

std::vector<double> differentiate_poly(const std::vector<double> & poly, int n)
{
  auto res = std::vector<double>(poly);
  int cs = poly.size();
  for (int i = 0; i < n; i++) {
    if ((cs - 1) == 0) {
      return {0.0};
    }
    
    cs -= 1;
    for (int i=0; i<cs; i++) {
      res[i] = res[i] * (double) (cs - i);
    }
  }

  res.resize(cs);

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
  for (size_t i = 1; i < in.size(); i++) {
    c = gcd(c, in[i]);
  }
  return c;
}

std::vector<int> reduce_gcd(std::vector<int> in)
{
  std::vector<int> out(in);
  int f = common_gcd(out);
  for (size_t i=0; i<out.size(); i++) {
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

void qmult_vec(const double * const As, const double * const B, 
               size_t n, double * const Cs)
{
  double Bs[16];
  Bs[0] = B[0];
  Bs[4] = -B[1];
  Bs[8] = -B[2];
  Bs[12] = -B[3];

  Bs[1] = B[1];
  Bs[5] = B[0];
  Bs[9] = B[3];
  Bs[13] = -B[2];

  Bs[2] = B[2];
  Bs[6] = -B[3];
  Bs[10] = B[0];
  Bs[14] = B[1];

  Bs[3] = B[3];
  Bs[7] = B[2];
  Bs[11] = -B[1];
  Bs[15] = B[0];

  dgemm_("N", "N", 4, n, 4, 1.0, Bs, 4, As, 4, 0.0, Cs, 4);
}

bool isclose(double a, double b)
{
  return fabs(a - b) <= (ATOL + RTOL * fabs(b));
}

void rotate_matrix(int m, int n, const double * const A,
                  const double * const B, double * C)
{
  double * temp = new double[m*n];
  
  // A_mn
  // B_nn
  // temp_nm = B_nn * A.T_nm
  // C_mm = A_mn * temp_nm

  dgemm_("T", "N", m, n, n, 1.0, A, n, B, n, 0.0, temp, m);
  dgemm_("N", "N", m, m, n, 1.0, temp, m, A, n, 0.0, C, m);

  delete [] temp;
}

int fact(int n)
{
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

double factorial(int n)
{
  return (double) fact(n);
}

void eigenvalues_sym(const double * const s, double * values)
{
  double F[9];
  usym(s, F);
  
  int swork = 15;
  double work[15];
  int info = 0;

  dsyev_("N", "U", 3, F, 3, values, work, swork, info);

  if (info != 0) throw LinalgError("Eigenvalue calculation failed");
}

void eigenvectors_sym(const double * const s, double * vectors)
{
  usym(s, vectors);

  int swork = 15;
  double work[15];
  int info = 0;
  double values[3];

  dsyev_("V", "U", 3, vectors, 3, values, work, swork, info);

  if (info != 0) throw LinalgError("Eigenvector calculation failed");
}

double I1(const double * const s)
{
  return s[0] + s[1] + s[2];
}

double I2(const double * const s)
{
  double F[9];
  usym(s, F);
  double tr1 = F[0] + F[4] + F[8];

  double F2[9];
  mat_mat(3, 3, 3, F, F, F2);

  double tr2 = F2[0] + F2[4] + F2[8];

  return 0.5 * (tr1*tr1 - tr2);
}

} // namespace neml
