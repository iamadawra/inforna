

#include "hairpin_energy.h"




/******************************************************
identifies the energy boni for special tetra loops
******************************************************/

double tetra_loop_energy(char* str)
{
   double energy = 0.0;

   //              Seq              Energy
   // --------------------------------------
   if      (strcmp(str,"GGGGAC") == 0) energy =  -3.00 ;
   else if (strcmp(str,"GGUGAC") == 0) energy =  -3.00 ;
   else if (strcmp(str,"CGAAAG") == 0) energy =  -3.00 ;
   else if (strcmp(str,"GGAGAC") == 0) energy =  -3.00 ;
   else if (strcmp(str,"CGCAAG") == 0) energy =  -3.00 ;
   else if (strcmp(str,"GGAAAC") == 0) energy =  -3.00 ;
   else if (strcmp(str,"CGGAAG") == 0) energy =  -3.00 ;
   else if (strcmp(str,"CUUCGG") == 0) energy =  -3.00 ;
   else if (strcmp(str,"CGUGAG") == 0) energy =  -3.00 ;
   else if (strcmp(str,"CGAAGG") == 0) energy =  -2.50 ;
   else if (strcmp(str,"CUACGG") == 0) energy =  -2.50 ;
   else if (strcmp(str,"GGCAAC") == 0) energy =  -2.50 ;
   else if (strcmp(str,"CGCGAG") == 0) energy =  -2.50 ;
   else if (strcmp(str,"UGAGAG") == 0) energy =  -2.50 ;
   else if (strcmp(str,"CGAGAG") == 0) energy =  -2.00 ;
   else if (strcmp(str,"AGAAAU") == 0) energy =  -2.00 ;
   else if (strcmp(str,"CGUAAG") == 0) energy =  -2.00 ;
   else if (strcmp(str,"CUAACG") == 0) energy =  -2.00 ;
   else if (strcmp(str,"UGAAAG") == 0) energy =  -2.00 ;
   else if (strcmp(str,"GGAAGC") == 0) energy =  -1.50 ;
   else if (strcmp(str,"GGGAAC") == 0) energy =  -1.50 ;
   else if (strcmp(str,"UGAAAA") == 0) energy =  -1.50 ;
   else if (strcmp(str,"AGCAAU") == 0) energy =  -1.50 ;
   else if (strcmp(str,"AGUAAU") == 0) energy =  -1.50 ;
   else if (strcmp(str,"CGGGAG") == 0) energy =  -1.50 ;
   else if (strcmp(str,"AGUGAU") == 0) energy =  -1.50 ;
   else if (strcmp(str,"GGCGAC") == 0) energy =  -1.50 ;
   else if (strcmp(str,"GGGAGC") == 0) energy =  -1.50 ;
   else if (strcmp(str,"GUGAAC") == 0) energy =  -1.50 ;
   else if (strcmp(str,"UGGAAA") == 0) energy =  -1.50 ;

   return energy;
};



/******************************************************************

finds best boni for tetraloops depending on the closing BP + the
first mismatch
     ==> fct. is not longer used when dealing with seq. constraints
*******************************************************************/
double best_tetra_loop_energy(int bp_i, int bp_j, int i, int j)
{
   double best_energy = 0.0;

   if ((bp_i == 0) && (bp_j == 3))
   {
      if ((i==2) && (j==0)) best_energy = -2.00;
   }
   else if ((bp_i == 1) && (bp_j == 2))
   {
      if ((i==2) && (j==0)) best_energy = -3.00;
      else if ((i==2) && (j==2)) best_energy = -2.5;
      else if ((i==3) && (j==1)) best_energy = -2.0;
      else if ((i==3) && (j==2)) best_energy = -3.0;
   }
   else if ((bp_i == 2) && (bp_j == 1))
   {
      if ((i==2) && (j==0)) best_energy = -3.00;
      else if ((i==2) && (j==2)) best_energy = -1.5;
      else if ((i==3) && (j==0)) best_energy = -1.5;
   }
   else if ((bp_i == 3) && (bp_j == 0))
   {
      if ((i==2) && (j==0)) best_energy = -1.50;
   }
   else if ((bp_i == 3) && (bp_j == 2))
   {
      if ((i==2) && (j==0)) best_energy = -2.50;
   }

   return best_energy;
};

/*************************************************************

calculation of the minimal free energy for an HL with given
size and given closing BP

*************************************************************/

double BestHairpinLoopEnergy(int bp_pos, int size, int bp_i, int bp_j)
{
   double energy = 0.0;
   double min = MAX_DOUBLE;
   double energy_help;

   int bp_pos_i = BP_Order[bp_pos][0];
   int bp_pos_j = BP_Order[bp_pos][1];

   // loop_destabilizing_energies
   //****************************
   if (size <= 30)
      energy += loop_destabilizing_energies[3*size-1];
   else
   {
      energy += loop_destabilizing_energies[3*30-1];
      energy += 1.75*RT*log((double)(size/30.0));
   }

   // terminal mismatch (size > 3)
   //**************************************
   if (size > 3)
   {
      // The exact assignment of the interior loop bases (not adjacent to the closing BP) is not important here.
      // Their valid assignments can be restricted by the constraints, but one assignment per position has to be
      // valid und this one is chosen for the best energy.

      // If the HPLoop would only consist of C's, a penalty of 0.3*size+1.6 has to be added.
      // This means that just in case that all bases in the loop are restricted to C by the
      // constraints the energy has to be corrected.
      // NOT TAKEN INTO ACCOUNT, SINCE THIS IS NOT TAKEN INTO ACCOUNT IN THE VIENNA PACKAGE AS WELL.

      // bool onlyCs = true;
      // for (int p=bp_pos_i+1; p<bp_pos_j; p++)
      //    if ( !((seq_constraints[p][1] == 1) && (Sum_SeqConst(p) == 1)))
      //    {
      //       onlyCs = false;
      //       break;
      //    }

      // if (onlyCs)
      //    energy = Sum_MaxDouble(energy,0.3*size+1.6);

      if (size != 4)
      {
         // finding the minimal terminal-mismatch energy for the given closing BP,
         // the remaining bases in the loop don't contribute to the energy, i.e. they
         // are not fixed until the traceback, their penalty is also considered there
         // (since at least one base has to be valid for each position, the exact
         // assignment has not to be fixed here)
         min = MAX_DOUBLE;
         for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
            {
               energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i), BasePenalty(bp_pos_j-1,j), mismatch_energies_hairpin[64*bp_i+16*i+4*bp_j+j]);
               if (energy_help < min)
                  min = energy_help;
            }
         energy = Sum_MaxDouble(energy,min);
      }
      else
      {
         char* tetra_plus_closing;
         tetra_plus_closing = (char*) malloc(sizeof(char)*6);
         tetra_plus_closing[0] = int2char(bp_i);
         tetra_plus_closing[5] = int2char(bp_j);

         // consider all cases of the tetraloop (4*4*4*4) and add term-mismatch and
         // possibly a bonus, furthermore add the penalties for all 4 bases
         //      i2   j2
         //   i           j
         //     bp_i-bp_j
         //       .   .

         min = MAX_DOUBLE;
         for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
               for (int i2=0; i2<4; i2++)
                  for (int j2=0; j2<4; j2++)
                  {
                     energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,i), BasePenalty(bp_pos_i+2,i2), BasePenalty(bp_pos_j-2,j2),BasePenalty(bp_pos_j-1,j));
                     tetra_plus_closing[1] = int2char(i);
                     tetra_plus_closing[2] = int2char(i2);
                     tetra_plus_closing[3] = int2char(j2);
                     tetra_plus_closing[4] = int2char(j);
                     energy_help = Sum_MaxDouble3(energy_help, tetra_loop_energy(tetra_plus_closing), mismatch_energies_hairpin[64*bp_i+16*i+4*bp_j+j]);

                     if (energy_help < min)
                        min = energy_help;
                  }
         energy = Sum_MaxDouble(energy,min);
         free(tetra_plus_closing);
      }
   }
   else if (size == 3)
   {
      min = MAX_DOUBLE;
      for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
            for (int i2=0; i2<4; i2++)
            {
               energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i), BasePenalty(bp_pos_i+2,i2), BasePenalty(bp_pos_j-1,j));

               // If the HL would only consist of Cs, a penalty of 1.4 has to be added.
               // This is just the best solution, if all other bases are forbidden by the constraints.
               if ((i==1) && (i2==1) && (j==1))
                  energy_help = Sum_MaxDouble(energy_help,Ctriloop);

               // If the HL would only consist of Gs, a term of -2.2 has to be added.
               // This is not considered in the ViennaPackage, so do we
               //if ((i==2) && (i2==2) && (j==2))
               //   energy_help = Sum_MaxDouble(energy_help,Gtriloop);

               if (energy_help < min)
                  min = energy_help;
            }
      energy = Sum_MaxDouble(energy,min);

      if (((bp_i==0) && (bp_j==3)) || ((bp_i==3) && (bp_j==0))|| ((bp_i==2) && (bp_j==3)) || ((bp_i==3) && (bp_j==2)))
         energy = Sum_MaxDouble(energy,terminalAU);
   }
   else
   {
      cerr << "HairpinLoopSize too small" << endl;
      exit(1);
   }
   return energy;
}




/*************************************************************

calculates energy fractions of a HL that are important for 
the difference of two assignments of it
(size and closing BP are given)

*************************************************************/

double HairpinLoopEnergy(int size, int pos_i, int* int_seq)
{
   double energy = 0.0;
   int i, k; // only_Cs;
   char* tetra_loop_assign;
   int pos_j = pos_i + size + 1;

   int bp_i = int_seq[pos_i];
   int bp_j = int_seq[pos_j];
   int mis_i = int_seq[pos_i+1];
   int mis_j = int_seq[pos_j-1];

   double energy_help = 0.0;

   // loop_destabilizing_energies (not necessary here, since this function is
   // just used for the difference of two assignments of the same HL)
   //------------------------------
   /*if (size <= 30)
      energy += loop_destabilizing_energies[3*size-1];
   else
   {
      energy += loop_destabilizing_energies[3*30-1];
      energy += 1.75*RT*log((double)(size/30.0));
   }*/


   // terminal mismatch (fuer size > 3)
   //----------------------------------
   if (size > 3)
   {
      //**************************************************************
      // If the HPLoop would only consist of C's, a penalty of 0.3*size+1.6 has to be added.
      // NOT TAKEN INTO ACCOUNT, SINCE THIS IS NOT TAKEN INTO ACCOUNT IN THE VIENNA PACKAGE AS WELL.
      /*only_Cs = 1;
      for (i=1; i<=size; i++)
         if (int_seq[pos_i+i] != 1)
         {
            only_Cs = 0;
            break;
         }
      if ((only_Cs == 1) && (size > 3))
         energy += 0.3*size+1.6;*/
      //**************************************************************

      //penalties for all bases of the loop
      energy_help = 0.0;
      for (i=1; i<=size; i++)
         energy_help = Sum_MaxDouble(energy_help,BasePenalty(pos_i+i,int_seq[pos_i+i]));

      // terminal-mismatch energy for the closing BP and its penalty
      energy = Sum_MaxDouble4(energy, energy_help, mismatch_energies_hairpin[64*bp_i+16*mis_i+4*bp_j+mis_j], PairPenalty(pos_i,pos_j,bp_i,bp_j));

      if (size == 4)
      {
         /*for tetraloops, consider a special term*/
         tetra_loop_assign = (char*) malloc(sizeof(char)*7);
         for (k=0; k<6; k++)
            tetra_loop_assign[k] = int2char(int_seq[pos_i+k]);
         tetra_loop_assign[6] = '\0';
         energy = Sum_MaxDouble(energy,tetra_loop_energy(tetra_loop_assign));
      }
   }
   else if (size == 3)
   {
      energy_help = Sum_MaxDouble3(BasePenalty(pos_i+1,mis_i), BasePenalty(pos_i+2,int_seq[pos_i+2]), BasePenalty(pos_j-1,mis_j));
      energy = Sum_MaxDouble3(energy, energy_help, PairPenalty(pos_i,pos_j,bp_i,bp_j));

      if (((bp_i==0) && (bp_j==3)) || ((bp_i==3) && (bp_j==0)) || ((bp_i==2) && (bp_j==3)) || ((bp_i==3) && (bp_j==2)))
         energy = Sum_MaxDouble(energy,terminalAU);
      // If the HL would only consist of Cs, a penalty of 1.4 has to be added.
      if ((int_seq[pos_i+1] == 1) && (int_seq[pos_i+2] == 1) && (int_seq[pos_i+3] == 1))
         energy = Sum_MaxDouble(energy,Ctriloop);
      // If the HL would only consist of Gs, a term of -2.2 has to be added. This is not considered in the ViennaPackage, so do we
      //if ((int_seq[pos_i+1] == 2) && (int_seq[pos_i+2] == 2) && (int_seq[pos_i+3] == 2))
      //   energy = Sum_MaxDouble(energy,Gtriloop);
   }
   else
   {
      printf("HairpinLoopSize too small\n");
      exit(1);
   }

   //if closingHL the last BP in a stack
   energy = Sum_MaxDouble(energy, Zero_or_StemEndAU(BP_Pos_Nr[pos_i], BP2int(bp_i,bp_j)));

   return energy;
}

