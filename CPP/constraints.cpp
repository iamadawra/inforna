
#include "constraints.h"

/*********************************************************
 tests, whether a sequence has a valid IUPAC-Codes,
 make T to U!!
*********************************************************/

int Check_iu()
{
   for(int i=0; i<(int)strlen(iupac_const); i++)
   {
      iupac_const[i] = toupper(iupac_const[i]);
      if (iupac_const[i] == 'T')
         iupac_const[i] = 'U';
      if  (Valid_IUPAC(iupac_const[i]) == 0)
         return 0;
   }
   return 1;
}



/*********************************************************
 tests, whether sequence constraints are possible to be
 fulfilled
*********************************************************/

int Check_constraints_bp()
{
   int bp_pos_i, bp_pos_j, bp_assign_i, bp_assign_j;
   bool possible;

   for (int bp_pos=0; bp_pos<numBP; bp_pos++)
   {
      bp_pos_i = BP_Order[bp_pos][0];
      bp_pos_j = BP_Order[bp_pos][1];

      possible = false;
      for (int bp_assign = 0; bp_assign<6; bp_assign++)
      {
         BP2_2(bp_assign, bp_assign_i, bp_assign_j);
         if (Compare_IUPAC_Base(iupac_const[bp_pos_i],bp_assign_i)+Compare_IUPAC_Base(iupac_const[bp_pos_j],bp_assign_j) == 2)
         {
            possible = true;
            break;
         }
      }
      if (possible == false)
      {
         cerr << "\n" << iupac_const[bp_pos_i] << " and " << iupac_const[bp_pos_j] << " are not compatible!\n\n";
         exit(1);
      }
   }
   return 1;
}



/***********************************************************
translates the IUPAC-Seq (iupac_const) the two dimen. array
where the constraints are stored, e.g. R becomes: 1010, i.e. 
A and G are valid, but C and U not
***********************************************************/

void getSeqConstraints()
{
   //allocate and set
   seq_constraints = (int**) malloc(sizeof(int*)*strlen(iupac_const));
   for (int i=0; i<(int)strlen(iupac_const); i++)
   {
      seq_constraints[i] = (int*) malloc(sizeof(int)*4);
      for (int j=0; j<4; j++)
         seq_constraints[i][j] = Compare_IUPAC_Base(iupac_const[i], j);
   }

   //print
   /*printf("seq_constraints: \n");
   for (int i=0; i<(int)strlen(iupac_const); i++)
   {
      for (int j=0; j<4; j++)
         printf("%d ",seq_constraints[i][j]);
      printf("\n");
   }*/
}

/***********************************************************
 tests, whether one of the values is MIN_INT
 if so, the sum is set to MIN_INT
***********************************************************/

double Sum_MinInt(double sum1, double sum2)
{
   if (((int)sum1 == MIN_INT) || ((int)sum2 == MIN_INT))
      return MIN_INT;
   else
      return sum1+sum2;
}


/***********************************************************
 tests, whether one of two values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

double Sum_MaxDouble(double sum1, double sum2)
{
   if ((sum1 == MAX_DOUBLE) || (sum2 == MAX_DOUBLE))
      return MAX_DOUBLE;
   else
      return sum1+sum2;
}


/***********************************************************
 tests, whether one of three values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

double Sum_MaxDouble3(double sum1, double sum2, double sum3)
{
   if ((sum1 == MAX_DOUBLE) || (sum2 == MAX_DOUBLE) || (sum3 == MAX_DOUBLE))
      return MAX_DOUBLE;
   else
      return sum1+sum2+sum3;
}


/***********************************************************
 tests, whether one of four values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

double Sum_MaxDouble4(double sum1, double sum2, double sum3, double sum4)
{
   if ((sum1 == MAX_DOUBLE) || (sum2 == MAX_DOUBLE) || (sum3 == MAX_DOUBLE) || (sum4==MAX_DOUBLE))
      return MAX_DOUBLE;
   else
      return sum1+sum2+sum3+sum4;
}


/***********************************************************
 tests, whether one of the values is MIN_INT
 if so, the difference is set to MIN_INT as well
 (used during the calculation of Ediff, is a constellation
 is penalized, this should be seen in Ediff as well)
***********************************************************/

double Sub_MinInt(double sum, double sub)
{
   if ((sum == MAX_DOUBLE) || (sub == MAX_DOUBLE))
      return MIN_INT;
   else
      return sum-sub;
}


/***********************************************************
 identifies, whether a given base is valid at a given pos.
 (if necessary, the penalty is given)
***********************************************************/

double BasePenalty(int pos, int base_assign)
{
   if (step == 1)
   {
      if (seq_constraints[pos][base_assign] == 1)
         return 0;
      else
         return MAX_DOUBLE;
   }
   else
   {
      if (seq_constraints[pos][base_assign] == 1)
         return 0;
      else if ((seq_constraints[pos][base_assign] == 0) && (mis_vec[pos] == 1))
         return 0;
      else
         return MAX_DOUBLE;
   }
}


/************************************************************
 identifies the penalty for a base pair
************************************************************/

double PairPenalty(int bp_pos, int bp_i, int bp_j)
{
   int pos_i = BP_Order[bp_pos][0];
   int pos_j = BP_Order[bp_pos][1];
   return Sum_MaxDouble(BasePenalty(pos_i, bp_i),BasePenalty(pos_j, bp_j));
}

double PairPenalty(int pos_i, int pos_j, int bp_i, int bp_j)
{
   return Sum_MaxDouble(BasePenalty(pos_i, bp_i),BasePenalty(pos_j, bp_j));
}

/*************************************************************
 sums up the sequence constraints for one position
 (= gives the number of valid bases)
*************************************************************/

int Sum_SeqConst(int pos_row)
{
   return seq_constraints[pos_row][0]+seq_constraints[pos_row][1]+seq_constraints[pos_row][2]+seq_constraints[pos_row][3];
}



/*************************************************************
 identifies a random free base
 (minds the constraints and free_bases_set2X)
*************************************************************/

int SetFreeBase(int pos)
{
   // if all free bases that give no energy-part are choosen to be set to a fixed base, we have to test the constraints at this positions 
   // nevertheless. If the fixed base is forbidden by the constraints, a base is chosen ramdomly among the other bases

   // generate a vector that gives the base assignments that are allowed (concerning free_bases_set2X and seq_constraints)
   int* free_bases_set2;
   free_bases_set2 = (int*) malloc(sizeof(int)*4);
   free_bases_set2[0] = Minimum(free_bases_set2A,seq_constraints[pos][0]);
   free_bases_set2[1] = Minimum(free_bases_set2C,seq_constraints[pos][1]);
   free_bases_set2[2] = Minimum(free_bases_set2G,seq_constraints[pos][2]);
   free_bases_set2[3] = Minimum(free_bases_set2U,seq_constraints[pos][3]);

   // If all bases fixed by free_bases_set2X are forbidden by seq_constraints, or if no fixed forces are given (free_bases_set2X = 0 for all X), 
   // this vector consists only of zeros and the "else" is done
   int sum_constraints = SumVec(free_bases_set2,4);

   if (sum_constraints > 0)
   {
      int rand_free = RandomBase(sum_constraints) + 1;
      int ones = 0;
      int column = -1;
      while (ones < rand_free)
      {
         column++;
         if (free_bases_set2[column] == 1)
            ones++;
      }
      //column is the randomly chosen base assignment (0=A, 1=C, 2=G, 3=U)
      return column;
   }
   else // if no forces because of free_bases_set2X or if all bases fixed by free_bases_set2X are forbidden by seq_constraints
   {
      int sum = Sum_SeqConst(pos); //number of valid assignments
      if (sum == 0)
      {
         cerr << "No valid base at position " << pos << "!\n";
         exit(1);
      }
      int rand = RandomBase(sum) + 1;  // the rand-th valid base assignment is chosen
      int ones = 0;
      int column = -1;
      while (ones < rand)
      {
         column++;
         if (seq_constraints[pos][column] == 1)
            ones++;
      }
      //column is the randomly chosen base assignment (0=A, 1=C, 2=G, 3=U)
      return column;
   }
}


