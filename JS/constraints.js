/*********************************************************
 tests, whether a sequence has a valid IUPAC-Codes,
 make T to U!!
*********************************************************/

function Check_iu()
{
   for(var i=0; i<iupac_const.length; i++)
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

function Check_constraints_bp()
{
   var bp_pos_i, bp_pos_j, bp_assign_i, bp_assign_j;
   bool possible;

   for (var bp_pos=0; bp_pos<numBP; bp_pos++)
   {
      bp_pos_i = BP_Order[bp_pos][0];
      bp_pos_j = BP_Order[bp_pos][1];

      possible = false;
      for (var bp_assign = 0; bp_assign<6; bp_assign++)
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
         throw("\n" + iupac_const[bp_pos_i] + " and " + iupac_const[bp_pos_j] + " are not compatible!\n\n");
      }
   }
   return 1;
}



/***********************************************************
translates the IUPAC-Seq (iupac_const) the two dimen. array
where the constraints are stored, e.g. R becomes: 1010, i.e. 
A and G are valid, but C and U not
***********************************************************/

function getSeqConstraints()
{
   //allocate and set
   seq_constraints = new Array((iupac_const).length);
   for (var i=0; i<(iupac_const.length); i++)
   {
      seq_constraints[i] = new Array(4);
      for (var j=0; j<4; j++)
         seq_constraints[i][j] = Compare_IUPAC_Base(iupac_const[i], j);
   }  
}

/***********************************************************
 tests, whether one of the values is MIN_INT
 if so, the sum is set to MIN_INT
***********************************************************/

function Sum_MinInt(sum1, sum2)
{
   if ((sum1 == -1000000) || ((int)sum2 == -1000000))
      return -1000000;
   else
      return sum1+sum2;
}


/***********************************************************
 tests, whether one of two values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

function Sum_MaxDouble(sum1, sum2)
{
   if ((sum1 == 1000000.0) || (sum2 == 1000000.0))
      return 1000000.0;
   else
      return sum1+sum2;
}


/***********************************************************
 tests, whether one of three values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

function Sum_MaxDouble3(sum1, sum2, sum3)
{
   if ((sum1 == 1000000.0) || (sum2 == 1000000.0) || (sum3 == 1000000.0))
      return 1000000.0;
   else
      return sum1+sum2+sum3;
}


/***********************************************************
 tests, whether one of four values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

function Sum_MaxDouble4(sum1, sum2, sum3, sum4)
{
   if ((sum1 == 1000000.0) || (sum2 == 1000000.0) || (sum3 == 1000000.0) || (sum4==1000000.0))
      return 1000000.0;
   else
      return sum1+sum2+sum3+sum4;
}


/***********************************************************
 tests, whether one of the values is MIN_INT
 if so, the difference is set to MIN_INT as well
 (used during the calculation of Ediff, is a constellation
 is penalized, this should be seen in Ediff as well)
***********************************************************/

function Sub_MinInt(sum, sub)
{
   if ((sum == 1000000.0) || (sub == 1000000.0))
      return -1000000;
   else
      return sum-sub;
}


/***********************************************************
 identifies, whether a given base is valid at a given pos.
 (if necessary, the penalty is given)
***********************************************************/

function BasePenalty(pos, base_assign)
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
         return 1000000.0;
   }
}


/************************************************************
 identifies the penalty for a base pair
************************************************************/

function PairPenalty(bp_pos, bp_i, bp_j)
{
   int pos_i = BP_Order[bp_pos][0];
   int pos_j = BP_Order[bp_pos][1];
   return Sum_MaxDouble(BasePenalty(pos_i, bp_i),BasePenalty(pos_j, bp_j));
}

function PairPenalty(pos_i, pos_j, bp_i, bp_j)
{
   return Sum_MaxDouble(BasePenalty(pos_i, bp_i),BasePenalty(pos_j, bp_j));
}

/*************************************************************
 sums up the sequence constraints for one position
 (= gives the number of valid bases)
*************************************************************/

function Sum_SeqConst(pos_row)
{
   return seq_constraints[pos_row][0]+seq_constraints[pos_row][1]+seq_constraints[pos_row][2]+seq_constraints[pos_row][3];
}



/*************************************************************
 identifies a random free base
 (minds the constraints and free_bases_set2X)
*************************************************************/

function SetFreeBase(pos)
{
   // if all free bases that give no energy-part are choosen to be set to a fixed base, we have to test the constraints at this positions 
   // nevertheless. If the fixed base is forbidden by the constraints, a base is chosen ramdomly among the other bases

   // generate a vector that gives the base assignments that are allowed (concerning free_bases_set2X and seq_constraints)
   var free_bases_set2;
   free_bases_set2 = new Array(4);
   free_bases_set2[0] = Minimum(free_bases_set2A,seq_constraints[pos][0]);
   free_bases_set2[1] = Minimum(free_bases_set2C,seq_constraints[pos][1]);
   free_bases_set2[2] = Minimum(free_bases_set2G,seq_constraints[pos][2]);
   free_bases_set2[3] = Minimum(free_bases_set2U,seq_constraints[pos][3]);

   // If all bases fixed by free_bases_set2X are forbidden by seq_constraints, or if no fixed forces are given (free_bases_set2X = 0 for all X), 
   // this vector consists only of zeros and the "else" is done
   var sum_constraints = SumVec(free_bases_set2,4);

   if (sum_constraints > 0)
   {
      var rand_free = RandomBase(sum_constraints) + 1;
      var ones = 0;
      var column = -1;
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
      var sum = Sum_SeqConst(pos); //number of valid assignments
      if (sum == 0)
      {
         throw("No valid base at position " + pos + "!\n");
      }
      var rand = RandomBase(sum) + 1;  // the rand-th valid base assignment is chosen
      var ones = 0;
      var column = -1;
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


