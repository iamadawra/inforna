#ifndef _INTERIOR_ENERGY__
#define _INTERIOR_ENERGY__

#include <stdlib.h>
#include "basics.h"
#include "constraints.h"

using namespace std;

double BestInteriorLoopEnergy(int bp_pos, int leftSize, int rightSize, int bp_i, int bp_j, int bp_before);
double InteriorLoopEnergy(int leftSize, int rightSize, int* int_seq, int pos_i, int pos_j);

#endif   // _INTERIOR_ENERGY_
