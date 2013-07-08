

#include "interior_energy.h"

/************************************************************

finds the minimal energy of an interior loop, if the
assignments of both closing BPs are known

************************************************************/

double BestInteriorLoopEnergy(int bp_pos, int leftSize, int rightSize, int bp_i, int bp_j, int bp_before)
{
   double energy = 0.0;
   double asym = 0.0;
   double min = MAX_DOUBLE;
   double energy_help;

   int size = leftSize + rightSize;
   int bp_new = BP2int(bp_i, bp_j);
   int base_i_before, base_j_before;
   int bp_pos_i, bp_pos_j;

   BP2_2(bp_before,base_i_before,base_j_before);
   bp_pos_i = BP_Order[bp_pos][0];
   bp_pos_j = BP_Order[bp_pos][1];

   // no bulges and loops of size 0:
   //******************************************
   if ((leftSize == 0) || (rightSize == 0))
   {
      cerr << "No interiorLoop!" << endl;
      exit(1);
   }

   // special cases:
   //******************************************
   if ((leftSize == 1) && (rightSize == 1))
   {
      for (int x=0; x<4; x++)
         for (int y=0; y<4; y++)
         {
            energy_help = Sum_MaxDouble3(interior_loop_1_1_energy[96*bp_new+24*x+4*bp_before+y], BasePenalty(bp_pos_i+1,x), BasePenalty(bp_pos_j-1,y));
            if (energy_help < min)
               min = energy_help;
         }
      energy = min;
   }
   else if ((leftSize == 1) && (rightSize == 2)) /*x(i-1)-x(i)=2 and y(i)-y(i-1)=3*/
   {
      for (int x=0; x<4; x++)
         for (int y=0; y<4; y++)
            for (int z=0; z<4; z++)
            {
               energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,x), BasePenalty(bp_pos_j-1,y), BasePenalty(bp_pos_j-2,z));
               energy_help = Sum_MaxDouble(energy_help, interior_loop_1_2_energy[384*bp_new+96*z+24*x+4*bp_before+y]);
               if (energy_help < min)
                  min = energy_help;
            }
      energy = min;
   }
   else if ((leftSize == 2) && (rightSize == 1)) /*x(i-1)-x(i)=3 and y(i)-y(i-1)=2*/
   {
      for (int x=0; x<4; x++)
         for (int y=0; y<4; y++)
            for (int z=0; z<4; z++)
            {
               energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,z), BasePenalty(bp_pos_i+2,y), BasePenalty(bp_pos_j-1,x));
               energy_help = Sum_MaxDouble(energy_help, interior_loop_1_2_energy[384*bp_before+96*z+24*x+4*bp_new+y]);
               if (energy_help < min)
                  min = energy_help;
            }
      energy = min;
   }
   else if ((leftSize == 2) && (rightSize == 2))
   {
      for (int x1=0; x1<4; x1++)
         for (int x2=0; x2<4; x2++)
            for (int y1=0; y1<4; y1++)
               for (int y2=0; y2<4; y2++)
               {
                  energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,x1), BasePenalty(bp_pos_i+2,y1), BasePenalty(bp_pos_j-1,x2), BasePenalty(bp_pos_j-2,y2));
                  energy_help = Sum_MaxDouble(energy_help, interior_loop_2_2_energy[1536*bp_new+256*bp_before+64*x1+16*x2+4*y1+y2]);
                  if (energy_help < min)
                     min = energy_help;
               }
      energy = min;
   }
   else
   {
      // size:
      //**********
      if (size <= 30)
         energy += loop_destabilizing_energies[3*(size-1)];
      else
      {
         energy += loop_destabilizing_energies[3*(30-1)];
         energy += 1.75*RT*log((double)(size/30.0));
      }

      // terminal mismatches a both closings:
      //***************************************
      // Attention: there are dependencies, if the IL has size 1 at one side of the loop,
      //            then this base is involved in the terminal mismatches at both closings
      if (leftSize == 1)      //        base_i_before - base_j_before
      {                       //  i1                                      i2
                              //                                          i3
                              //                 bp_i - bp_j
         for (int i1=0; i1<4; i1++)
            for (int i2=0; i2<4; i2++)
               for (int i3=0; i3<4; i3++)
               {
                  energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i1), BasePenalty(bp_pos_j-rightSize,i2), BasePenalty(bp_pos_j-1,i3));
                  energy_help = Sum_MaxDouble3(energy_help, mismatch_energies_interior[64*base_j_before+16*i2+4*base_i_before+i1], mismatch_energies_interior[64*bp_i+16*i1+4*bp_j+i3]);

                  if (energy_help < min)
                     min = energy_help;
               }
         energy = Sum_MaxDouble(energy, min);
      }
      else if (rightSize == 1)//        base_i_before - base_j_before
      {                       //  i1                                      i3
                              //  i2
                              //                 bp_i - bp_j
         for (int i1=0; i1<4; i1++)
            for (int i2=0; i2<4; i2++)
               for (int i3=0; i3<4; i3++)
               {
                  energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+leftSize,i1), BasePenalty(bp_pos_i+1,i2), BasePenalty(bp_pos_j-1,i3));
                  energy_help = Sum_MaxDouble3(energy_help, mismatch_energies_interior[64*base_j_before+16*i3+4*base_i_before+i1], mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i3]);

                  if (energy_help < min)
                     min = energy_help;
               }
         energy = Sum_MaxDouble(energy,min);
      }
      else                    //        base_i_before - base_j_before
      {                       //  i1                                      i3
                              //  i2                                      i4
                              //                 bp_i - bp_j
           for (int i1=0; i1<4; i1++)
              for ( int i3=0; i3<4; i3++)
              {
                 energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+leftSize,i1), BasePenalty(bp_pos_j-rightSize,i3), mismatch_energies_interior[64*base_j_before+16*i3+4*base_i_before+i1]);
                 if (energy_help < min)
                    min = energy_help;
              }
           energy = Sum_MaxDouble(energy,min);

           min = MAX_DOUBLE;
           for (int i2=0; i2<4; i2++)
              for ( int i4=0; i4<4; i4++)
              {
                 energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i2), BasePenalty(bp_pos_i-1,i4), mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i4]);
                 if (energy_help < min)
                    min = energy_help;
              }
           energy = Sum_MaxDouble(energy,min);
      }

      // asymmetry-penalty:
      //********************
      // (for 1x2 IL not extra penalty)
      if (leftSize != rightSize)
      {
         asym = 0.5 * abs(leftSize-rightSize);
         if (asym > 3.0)
            asym = 3.0;
         energy = Sum_MaxDouble(energy,asym);
      }
   }
   return energy;
}



/************************************************************

finds the energy of an interior loop, if the assignments of
both closing BPs and all free bases are known

(only energy fragments that are important for the difference 
of two assignments of the loop)

************************************************************/

double InteriorLoopEnergy(int leftSize, int rightSize, int* int_seq, int pos_i, int pos_j)
{
   double energy = 0.0;
   //double asym = 0.0;
   int x1, x2, y1, y2, i1, i2, i3, i4;
   //int size = leftSize + rightSize;

   int bp_i = int_seq[pos_i];
   int bp_j = int_seq[pos_j];
   int bp_i_before = int_seq[pos_i+leftSize+1];
   int bp_j_before = int_seq[pos_j-rightSize-1];

   int bp_new = BP2int(bp_i, bp_j);
   int bp_before = BP2int(bp_i_before, bp_j_before);

   // no bulges and loops of size 0:
   //--------------------------------------------
   if ((leftSize == 0) || (rightSize == 0))
   {
      printf("No interiorLoop!\n");
      exit(1);
   }

   // Special cases:
   //------------------------------------------------------
   if ((leftSize == 1) && (rightSize == 1))
   {
      x1 = int_seq[pos_i+1];
      y1 = int_seq[pos_j-1];
      energy += interior_loop_1_1_energy[96*bp_new+24*x1+4*bp_before+y1];
   }
   else if ((leftSize == 1) && (rightSize == 2)) /*x(i-1)-x(i)=2 and y(i)-y(i-1)=3*/
   {
      x1 = int_seq[pos_i+1];
      y1 = int_seq[pos_j-2];
      y2 = int_seq[pos_j-1];
      energy += interior_loop_1_2_energy[384*bp_new+96*y1+24*x1+4*bp_before+y2];
   }
   else if ((leftSize == 2) && (rightSize == 1)) /*x(i-1)-x(i)=3 and y(i)-y(i-1)=2*/
   {
      x1 = int_seq[pos_i+1];
      x2 = int_seq[pos_i+2];
      y1 = int_seq[pos_j-1];
      energy += interior_loop_1_2_energy[384*bp_before+96*x1+24*y1+4*bp_new+x2];
   }
   else if ((leftSize == 2) && (rightSize == 2))
   {
      x1 = int_seq[pos_i+1];
      x2 = int_seq[pos_i+2];
      y1 = int_seq[pos_j-2];
      y2 = int_seq[pos_j-1];
      energy += interior_loop_2_2_energy[1536*bp_new+256*bp_before+64*x1+16*y2+4*x2+y1];
   }
   else
   {
      // size: (is not important here, since this function is just used for the calculation of a difference
      // of two assignments of the same IL)
      //----------------------------------------
      /*if (size <= 30)
         energy += loop_destabilizing_energies[3*(size-1)];
      else
      {
         energy += loop_destabilizing_energies[3*(30-1)];
         energy += 1.75*RT*log((double)(size/30.0));
      }*/

      // terminal mismatches a both closings:
      //----------------------------------------
      // Attention: there are dependencies, if the IL has size 1 at one side of the loop,
      //            then this base is involved in the terminal mismatches at both closings
      if (leftSize == 1)      //        bp_i_before - bp_j_before
      {                       //  i1                                      i2
                              //                                          i3
                              //               bp_i - bp_j
         i1 = int_seq[pos_i+1];
         i2 = int_seq[pos_j-rightSize];
         i3 = int_seq[pos_j-1];
         energy += mismatch_energies_interior[64*bp_j_before+16*i2+4*bp_i_before+i1] + mismatch_energies_interior[64*bp_i+16*i1+4*bp_j+i3];
      }
      else if (rightSize == 1)//        bp_i_before - bp_j_before
      {                       //  i1                                      i3
                              //  i2
                              //               bp_i - bp_j
         i1 = int_seq[pos_i+leftSize];
         i2 = int_seq[pos_i+1];
         i3 = int_seq[pos_j-1];
         energy += mismatch_energies_interior[64*bp_j_before+16*i3+4*bp_i_before+i1] + mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i3];
      }
      else                    //        bp_i_before - bp_j_before
      {                       //  i1                                      i3
                              //  i2                                      i4
                              //               bp_i - bp_j
           i1 = int_seq[pos_i+leftSize];
           i2 = int_seq[pos_i+1];
           i3 = int_seq[pos_j-rightSize];
           i4 = int_seq[pos_j-1];
           energy += mismatch_energies_interior[64*bp_j_before+16*i3+4*bp_i_before+i1] + mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i4];
      }

      // asymmetry-penalty: (not important here)
      //------------------------
      // (for 1x2 IL not extra penalty)
      /*if (leftSize != rightSize)
      {
         asym = 0.5 * abs(leftSize-rightSize);
         if (asym > 3.0)
            asym = 3.0;
         energy+=asym;
      }*/
   }
   return energy;
}



