#ifndef _HAIRPIN_ENERGY__
#define _HAIRPIN_ENERGY__

#include <stdlib.h>
#include "basics.h"
#include "constraints.h"

using namespace std;

double tetra_loop_energy(char* str);
double best_tetra_loop_energy(int bp_i, int bp_j, int i, int j);
double BestHairpinLoopEnergy(int bp_pos, int size, int bp_i, int bp_j);
double HairpinLoopEnergy(int size, int pos_i, int* int_seq);

#endif   // _HAIRPIN_ENERGY_

