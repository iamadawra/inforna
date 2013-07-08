#include "energy.h"

/******************************************************
tests, whether the end of a stack (stem) is reached
******************************************************/

bool StackEnd(int bp_pos)
{
   bool end = false;
   // if there exists a next BP
   if (bp_pos+1 < numBP)
   {
      // if y(i+1) < x(i) or if (x(i+1),y(i+1)) closing BP of a ML
      if ((BP_Order[bp_pos+1][1] < BP_Order[bp_pos][0]) || (BP_Order[bp_pos+1][3] != 0))
         end = true;
   }
   else //if it is the last BP, it is the last BP in the stem as well
      end = true;

   return end;
}


/******************************************************
A-U penalty at the end of a stem,

HERE: test, whether the end of a stem is reached and if
the BP is assigned to A-U or U-A or G-U or U-G,
if so, then return: 0.5 otherwise: 0
******************************************************/

double Zero_or_StemEndAU(int bp_pos, int bp_assign)
{
   if ((StackEnd(bp_pos)) && ((bp_assign==0) || (bp_assign==3) || (bp_assign==4) || (bp_assign==5)))
      return terminalAU;
   else
      return 0.0;
}

