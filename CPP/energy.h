#ifndef _ENERGY__
#define _ENERGY__

#include <stdlib.h>
#include "basics.h"
#include "constraints.h"

#include "./data/stacking_energies.dat"
#include "./data/interior_loop_1_1_energies.dat"
#include "./data/interior_loop_1_2_energies.dat"
#include "./data/interior_loop_2_2_energies.dat"
#include "./data/loop_destabilizing_energies.dat"
#include "./data/single_base_stacking_energies.dat"
#include "./data/terminal_mismatch_hairpin.dat"
#include "./data/terminal_mismatch_interior.dat"

using namespace std;

const double terminalAU = 0.5;
const double Ctriloop = 1.4;
const double Gtriloop = -2.2;

bool StackEnd(int bp_pos);
double Zero_or_StemEndAU(int bp_pos, int bp_assign);

#endif   // _ENERGY_
