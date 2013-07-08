
#include "stacking_energy.h"

double StackingEnergy(int bp_i, int bp_j, int bp_before)
{
   double energy = 0.0;
   int base_i_before, base_j_before;

   BP2_2(bp_before,base_i_before,base_j_before);
   energy += stacking_energies[64*bp_i+16*base_i_before+4*bp_j+base_j_before];

   return energy;
}


