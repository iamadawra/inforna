#ifndef _CONSTRAINTS__
#define _CONSTRAINTS__

#include <stdlib.h>
#include "basics.h"

using namespace std;

int Check_iu();
int Check_constraints_bp();
void getSeqConstraints();

double Sum_MinInt(double sum1, double sum2);
double Sum_MaxDouble(double sum1, double sum2);
double Sum_MaxDouble3(double sum1, double sum2, double sum3);
double Sum_MaxDouble4(double sum1, double sum2, double sum3, double sum4);
double Sub_MinInt(double sum, double sub);

double BasePenalty(int pos, int base_assign);
double PairPenalty(int bp_pos, int bp_i, int bp_j);
double PairPenalty(int pos_i, int pos_j, int bp_i, int bp_j);
int Sum_SeqConst(int pos_row);
int SetFreeBase(int pos);

#endif   // _CONSTRAINTS_
