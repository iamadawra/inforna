var MAXALPHA = 20;
var TIME_OUT_TIME = 3600;
var MIN_INT = - 9007199254740992;
var MAX_INT = 9007199254740992;

function PosEnergy(pos, energy) {
	this.pos = pos;
	this.energy = energy;

}

var base = 4, npairs = 6;

var  BP_Precursors;   // precursor for each BP (ML,EL have more than one, this is not stored here, instead set to -1)
var   BP_Successors;  // successor for each BP
var Ediff;        // energy diffence that arises if a free base or a BP is changed
var av_Ediff;      // average energy difference if a free base or a BP is changed (average over all possible changes per pos.)
var max_Ediff;     // maximal energy difference if a free base or a BP is changed (max. over all possible changes per pos.)

var time_out = 0;      // if the maximal running time is exceeded: set to 1
var start_time;
var zw_time;


/*-------------------------------------------------------------------------*/
var fold_type;
var cost2;

/*---------------------------------------------------------------------------*/

/**********************************************************
creates a vector in which for each pos. in brackets the
pos. in BP_Order is stored, free bases are set to -1
***********************************************************/

function Pos2BP_Pos()
{
   BP_Pos_Nr = new Array(struct_len);
   var i, bp;
   for (i=0; i<struct_len; i++)
      BP_Pos_Nr[i] = -1;

   for (bp=0; bp<numBP; bp++)
   {
      BP_Pos_Nr[BP_Order[bp][0]] = bp;
      BP_Pos_Nr[BP_Order[bp][1]] = bp;
   }
}

/*---------------------------------------------------------------------------*/
/*********************************************************
 finds the predecessors of a BP
 (a closing pair of a ML has more than one predecessors,
 these are not needed here, therefore: just set to -1,

 but all predecessors of a ML-BP are important for
 finding the successors, thus an exact array including
 all predecessors is stored addidionally (prec) )
*********************************************************/

function GetPrecursors()
{
   var i,j,vorg;
   var all_vorg;
   var prec;
   prec = new Array(numBP+1);
   for (i=0; i<(numBP+1); i++)
   {
      //allocate as much memory as many predecessors the BP has (at least one pos.)
      vorg = Math.max(1, BP_Order[i][3]);
      prec[i] = new Array(vorg);
      for (j=0; j<vorg; j++)
         prec[i][j] = -1;
   }

   BP_Precursors = new Array(numBP + 1);
   for (i=0; i<(numBP+1); i++)
      BP_Precursors[i] = -1;

   for (i=0; i<numBP; i++)
   {
      // one or no predecessor
      if (BP_Order[i][3] == 0)
      {
         if (i>0)
            if ((BP_Order[i][0] < BP_Order[i-1][0]) && (BP_Order[i][1] > BP_Order[i-1][1]))
            {
               prec[i][0] = i-1;
               BP_Precursors[i] = i-1;
            }
      }
      else
      {
         all_vorg = FindStemEnds(i);
         for (j=0; j<BP_Order[i][3]; j++)
            prec[i][j] = all_vorg[j];
         //for closingML -1 in BP_Precursors, also for dangling ends
      }
   }

   //external loop (dangling ends) extra
   j=0;
   i=0;
   while (j < brackets.lengtH)
   {
      if (brackets[j] == '(')
      {
         prec[numBP][i++] = BP_Pos_Nr[j];
         j = BP_Order[BP_Pos_Nr[j]][1] + 1;
      }
      else
         j++;
   }

   return prec;
}

/*---------------------------------------------------------------------------*/

/*********************************************************
 finding all successsors of the BPs
*********************************************************/

function GetSuccessors(precs)
{
   var i,j;

   BP_Successors = new Array(numBP+1);
   for (i=0; i<(numBP+1); i++)
      BP_Successors[i] = -1;

   for (i=0; i<(numBP+1); i++)
      for (j=0; j<Math.max(1,BP_Order[i][3]); j++)  //i has BP_Order[i][3] many predecessors
         if (precs[i][j] != -1)
            BP_Successors[precs[i][j]] = i;
}

/*---------------------------------------------------------------------------*/
/******************************************************
 allocates memory for Ediff
******************************************************/

function alloc_Ediff()
{
   var i, beleg;
   Ediff = new Array(struct_len);
   av_Ediff = new Array(struct_len);
   max_Ediff = new Array(struct_len);

   for (i=0; i<struct_len; i++)
   {
      if (BP_Pos_Nr[i] == -1)
         beleg = 4;
      else if (i == BP_Order[BP_Pos_Nr[i]][1]) //closing brackets
         beleg = 6;
      else
         beleg = 0;

      Ediff[i] = new Array(beleg);
   }
   //init_Ediff();
}

/*---------------------------------------------------------------------------*/
/***********************************************************
 initializing the Ediff array
************************************************************/

function init_Ediff()
{
   var i, j, beleg;

   for (i=0; i<struct_len; i++)
   {
      av_Ediff[i] = MIN_INT;
      max_Ediff[i] = MIN_INT;

      if (BP_Pos_Nr[i] == -1)
         beleg = 4;
      else if (i == BP_Order[BP_Pos_Nr[i]][1]) //closing bracket
         beleg = 6;
      else
         beleg = 0;

      for (j=0; j<beleg; j++)
         Ediff[i][j] = MIN_INT;
   }
}

/*---------------------------------------------------------------------------*/
/***********************************************************
printing the Ediff array
************************************************************/

function print_Ediff()
{
   var i, j, beleg;

   outln("Ediff: ");
   for (i=0; i<struct_len; i++)
   {
      if (BP_Pos_Nr[i] == -1)
         beleg = 4;
      else if (i == BP_Order[BP_Pos_Nr[i]][1]) //closing bracket
         beleg = 6;
      else
         beleg = 0;

      out(i + " ");
      for (j=0; j<beleg; j++)
         out(Ediff[i][j] + " ");
      out("\n");
   }
}
/*---------------------------------------------------------------------------*/
/***********************************************************
printing the av_Ediff array
************************************************************/

function print_av_Ediff()
{
   var i;
   outln("av_Ediff:");
   for (i=0; i<struct_len; i++)
   {
      out(i + " ");
      outln(av_Ediff[i]);
   }
}
/*---------------------------------------------------------------------------*/
/***********************************************************
calculating the av_Ediff array
************************************************************/

function make_av_Ediff()
{
   var i, j, beleg;
   var av;

   for (i=0; i<struct_len; i++)
   {
      av = 0.0;
      if (BP_Pos_Nr[i] == -1)
         beleg = 4;
      else if (i == BP_Order[BP_Pos_Nr[i]][1]) //closing bracket
         beleg = 6;
      else
         beleg = 0;

      for (j=0; j<beleg; j++)
         //av = Sum_MinInt(av, Ediff[i][j]);
    av += Ediff[i][j];

      if ((beleg != 0) && (av != MIN_INT))
         av_Ediff[i] = av / beleg;
      else
         av_Ediff[i] = MIN_INT;
   }
}
/*---------------------------------------------------------------------------*/
/***********************************************************
printing the max_Ediff array
************************************************************/

function print_max_Ediff()
{
   var i;
   outln("max_Ediff:");
   for (i=0; i<struct_len; i++)
   {
      out(i + " ");
      out(max_Ediff[i] + "\n");
   }
}
/*---------------------------------------------------------------------------*/
/***********************************************************
calculation the max_Ediff array
************************************************************/

function make_max_Ediff()
{
   var i, beleg;
   var max;

   for (i=0; i<struct_len; i++)
   {
      if (BP_Pos_Nr[i] == -1)
         beleg = 4;
      else if (i == BP_Order[BP_Pos_Nr[i]][1]) //closing bracket
         beleg = 6;
      else
         beleg = 0;

      if (beleg != 0)
      {
         max = MaxiVec(Ediff[i], beleg);
         max_Ediff[i] = max[1];
      }
      else
         max_Ediff[i] = MIN_INT;
   }
}

/*---------------------------------------------------------------------------*/


/***********************************************************
 Positions are sorted according to the amount of Ediff
************************************************************/

function compare_Pos(a, b)
{
   if(a.energy > b.energy) return -1;
   if(a.energy < b.energy) return 1;
   return 0;
}


/*---------------------------------------------------------------------------*/
/***********************************************************
 fixes mut_pos_list (order of mutation)
************************************************************/

function make_mut_pos_list()
{
   var list = new Array(struct_len);

   //make an array of PosEnergy-structs from max_Ediff
   var av_E = new Array(struct_len);

   for (var i=0; i<struct_len; i++)
   {
      av_E[i] = PosEnergy(i, av_Ediff[i]);
   }

   av_E.sort(compare_Pos);

   for (var i=0; i<struct_len; i++)
      list[i] = av_E[i].pos;

   return list;
}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/***********************************************************************
 determines the kind of the structural component and its energy
***********************************************************************/

function getPartEnergy(bp_pos, bp_pos_before, int_seq)
{
   var k, pos_i, pos_j, pos_before_i, pos_before_j, assign_i, assign_j, assign_before_i, assign_before_j, assign_bp, assign_bp_before;
   var energy_help = 0.0;

   pos_i = BP_Order[bp_pos][0];
   pos_j = BP_Order[bp_pos][1];
   pos_before_i = BP_Order[bp_pos_before][0];
   pos_before_j = BP_Order[bp_pos_before][1];

   assign_i = int_seq[pos_i];
   assign_j = int_seq[pos_j];
   assign_before_i = int_seq[pos_before_i];
   assign_before_j = int_seq[pos_before_j];

   //***********************
   //Stack
   //***********************
   if ((pos_before_i - pos_i == 1) && (pos_j - pos_before_j == 1))
   {
      assign_bp = BP2int(assign_i, assign_j);
      assign_bp_before = BP2int(assign_before_i,assign_before_j);
      //base pair penalty only for the second closing base pair
      //return Sum_MaxDouble4(StackingEnergy(assign_i, assign_j, assign_bp_before), PairPenalty(pos_i,pos_j,assign_i,assign_j), PairPenalty(pos_before_i, pos_before_j, assign_before_i, assign_before_j), Zero_or_StemEndAU(bp_pos, assign_bp));
      return Sum_MaxDouble3(StackingEnergy(assign_i, assign_j, assign_bp_before), PairPenalty(pos_i,pos_j,assign_i,assign_j), Zero_or_StemEndAU(bp_pos, assign_bp));
   }

   //***********************
   //Bulge (right)
   //***********************
   else if ((pos_before_i-pos_i == 1) && (pos_j - pos_before_j > 1))
   {
      assign_bp = BP2int(assign_i, assign_j);
      assign_bp_before = BP2int(assign_before_i,assign_before_j);
      //Penalties for all pos. in the bulge
      energy_help = 0.0;
      for (k=1; k<pos_j-pos_before_j; k++)
         energy_help = Sum_MaxDouble(energy_help, BasePenalty(pos_before_j+k, int_seq[pos_before_j+k])); 
      //base pair penalties
      //energy_help = Sum_MaxDouble3(energy_help, PairPenalty(pos_i,pos_j,assign_i,assign_j), PairPenalty(pos_before_i, pos_before_j, assign_before_i, assign_before_j));
      //base pair penalty only for the second closing base pair
      energy_help = Sum_MaxDouble(energy_help, PairPenalty(pos_i,pos_j,assign_i,assign_j));
      return Sum_MaxDouble3(energy_help,BulgeEnergy(pos_j-pos_before_j,assign_i, assign_j, assign_bp_before), Zero_or_StemEndAU(bp_pos, assign_bp));
   }

   //***********************
   //Bulge (left)
   //***********************
   else if ((pos_before_i-pos_i > 1) && (pos_j - pos_before_j == 1))
   {
      assign_bp = BP2int(assign_i, assign_j);
      assign_bp_before = BP2int(assign_before_i,assign_before_j);
      //Penalties for all pos. in the bulge
      energy_help = 0.0;
      for (k=1; k<pos_before_i-pos_i; k++)
         energy_help = Sum_MaxDouble(energy_help, BasePenalty(pos_i+k, int_seq[pos_i+k])); 
      //base pair penalties
      //energy_help = Sum_MaxDouble3(energy_help, PairPenalty(pos_i,pos_j,assign_i,assign_j), PairPenalty(pos_before_i, pos_before_j, assign_before_i, assign_before_j));
      //base pair penalty only for the second closing base pair
      energy_help = Sum_MaxDouble(energy_help, PairPenalty(pos_i,pos_j,assign_i,assign_j));
      return Sum_MaxDouble3(energy_help, BulgeEnergy(pos_before_i-pos_i,assign_i,assign_j,assign_bp_before),  Zero_or_StemEndAU(bp_pos, assign_bp));
   }

   //***********************
   //ML or dangl end.
   //***********************
   else if (BP_Order[bp_pos][3] != 0)
   {
      //ML
      if (bp_pos != numBP)
         return MLEnergy(bp_pos,int_seq);
      //dangling end
      else
         return externEnergy(int_seq);
   }

   //***********************
   //IL
   //***********************
   else
   {
      assign_bp = BP2int(assign_i, assign_j);
      assign_bp_before = BP2int(assign_before_i,assign_before_j);
      //test all bases in the loop, add penalties
      energy_help = 0.0;
      for (k=1; k<pos_before_i-pos_i; k++)
         energy_help = Sum_MaxDouble(energy_help, BasePenalty(pos_i+k, int_seq[pos_i+k]));
      for (k=1; k<pos_j-pos_before_j; k++)
         energy_help = Sum_MaxDouble(energy_help, BasePenalty(pos_before_j+k, int_seq[pos_before_j+k])); 
      //base pair penalties
      //energy_help = Sum_MaxDouble3(energy_help, PairPenalty(pos_i,pos_j,assign_i,assign_j), PairPenalty(pos_before_i, pos_before_j, assign_before_i, assign_before_j));
      //base pair penalty only for the second closing base pair
      energy_help = Sum_MaxDouble(energy_help, PairPenalty(pos_i,pos_j,assign_i,assign_j));
      return Sum_MaxDouble3(energy_help, InteriorLoopEnergy(pos_before_i-pos_i-1,pos_j-pos_before_j-1,int_seq,pos_i,pos_j), Zero_or_StemEndAU(bp_pos, assign_bp));
   }
}



/*---------------------------------------------------------------------------*/
/*************************************************************
 calculates the part of energy that is influenced by a BP
 (even if in mfe-mode only parts of the structures are analyzed,
 the energy difference refers to the whole structure)
*************************************************************/
function get_BP_Energy(pos_i, int_seq)
{
   var energy = 0.0;

   var pos_BP;     //pos_i is the absolute pos. in the structure, pos_BP is the pos. of the BP in BP_Order (in pos_i is a BP)
   var pos_j;      //binding pos. of pos_i, if there is a BP in pos_i
   var bp_i, bp_j; // assignment of the BP
   var size;       // size of the loop
   var pos_i_vor, pos_j_vor, pos_i_nach, pos_j_nach; // pos. of the previous and the following BPs (if there are exactly one prec. and one succ.)
   var pos_BP_vor, pos_BP_nach;                      // pos. of the previous and the following BP in BP_Order
   var bp_i_vor, bp_j_vor, bp_i_nach, bp_j_nach;     // assignments of the previous and the following BP

   pos_j = BP_Order[BP_Pos_Nr[pos_i]][1];
   bp_i = int_seq[pos_i];
   bp_j = int_seq[pos_j];

   pos_BP = BP_Pos_Nr[pos_i];
   pos_BP_vor = BP_Precursors[pos_BP];
   if (pos_BP_vor != -1)
   {
      pos_i_vor = BP_Order[pos_BP_vor][0];
      pos_j_vor = BP_Order[pos_BP_vor][1];
      bp_i_vor = int_seq[pos_i_vor];
      bp_j_vor = int_seq[pos_j_vor];
   }
   pos_BP_nach = BP_Successors[pos_BP];
   if (pos_BP_nach != -1)
   {
      pos_i_nach = BP_Order[pos_BP_nach][0];
      pos_j_nach = BP_Order[pos_BP_nach][1];
      bp_i_nach = int_seq[pos_i_nach];
      bp_j_nach = int_seq[pos_j_nach];
   }

   //depending on the kind of the BP:
   //=====================================
   //CLOSING HL or ML or dangling end:

   if (pos_BP_vor == -1)  // no predecessor
   {
      //CLOSING HL:
      //---------------
      if (BP_Order[pos_BP][3] == 0) // no ML
      {
         size = BP_Order[pos_BP][1] - BP_Order[pos_BP][0] - 1;
         energy = Sum_MaxDouble(energy, HairpinLoopEnergy(size,pos_i,int_seq));

         //if pos_BP_nach == -1, pos_BP would already be the imaginary last BP ("closing"
         //the external loop) and no further structural element follows
         if (pos_BP_nach != -1)
            energy = Sum_MaxDouble(energy, getPartEnergy(pos_BP_nach, pos_BP, int_seq));
      }

      //CLOSING ML or dangling end:
      //-------------------------------
      else
      {
         if (pos_BP_nach != -1)
         {
            energy = Sum_MaxDouble(energy, MLEnergy(pos_BP,int_seq));
            energy = Sum_MaxDouble(energy, getPartEnergy(pos_BP_nach, pos_BP, int_seq));
         }
         else  //dangling has no predecessor and no successor
         {
            energy = Sum_MaxDouble(energy, externEnergy(int_seq));
         }
      }
   }
   //BL or IL or successor is dangling end or closing ML:
   //-------------------------------------------------------
   else
   {
      //has predecessor
      energy = Sum_MaxDouble(energy, getPartEnergy(pos_BP, pos_BP_vor, int_seq));
      //all that have a predecessor, have a successor as well
      energy = Sum_MaxDouble(energy, getPartEnergy(pos_BP_nach, pos_BP, int_seq));
   }

   return energy;
}

/*---------------------------------------------------------------------------*/
/*************************************************************
 calculates the part of energy that is influenced by a free base
 (even if in mfe-mode only parts of the structures are analyzed,
 the energy difference refers to the whole structure)
*************************************************************/
function get_BasePart_Energy(pos_i, int_seq)
{
   var bp_i;                    // assignment of the base
   var size;                    // loop size
   var pos_BP_vor, pos_BP_nach; // pos. of the BPs that enclose the free base (in BP_Order)
   var i, j;
   var group_i, group_j;        // left and right border of the group of free bases where the consideres free base (pos_i) is located

   // ATTENTION!!: PREDECESSOR AND SUCCESSOR BP DON'T HAVE THE SAME MEANING AS FOR BPs! HERE THE BPs ARE NOT 
   //              PREDECESSOR AND SUCCESSOR TO EACH OTHER AUTOMATICALLY!

   bp_i = int_seq[pos_i];

   //finding the group of free bases
   //+++++++++++++++++++++++++++++++++++++++++++
   i = pos_i;
   while ((0<=i) && (brackets[i] == '.'))
      i--;
   group_i = i+1;

   j = pos_i;
   while ((j<struct_len) && (brackets[j] == '.'))
      j++;
   group_j = j-1;

   // if the free base is not adjacent to a stem, it has no influence on the energy
   // thus, return 0.0 + penalty
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if ((pos_i == 0) && (brackets[pos_i+1] == '.'))
      return Sum_MaxDouble(0.0,BasePenalty(pos_i, int_seq[pos_i]));
   else if ((pos_i == struct_len-1) && (brackets[pos_i-1] == '.'))
      return Sum_MaxDouble(0.0,BasePenalty(pos_i, int_seq[pos_i]));
   else if ((brackets[pos_i-1] == '.') && (brackets[pos_i+1] == '.'))
      return Sum_MaxDouble(0.0,BasePenalty(pos_i, int_seq[pos_i]));

   else
   {
      // EL (- G i) or (i G -)
      if ((group_i == 0) || (group_j == struct_len-1))
      {
         return externEnergy(int_seq);
      }

      pos_BP_vor = BP_Pos_Nr[group_i-1];
      pos_BP_nach = BP_Pos_Nr[group_j+1];

      // in HairpinLoop
      if (pos_BP_vor == pos_BP_nach)
      {
         size = BP_Order[pos_BP_vor][1] - BP_Order[pos_BP_vor][0] - 1;
         return HairpinLoopEnergy(size,BP_Order[pos_BP_vor][0],int_seq);
      }

      // free base is right in a BL or IL or ML, pos_BP_nach = closingML
      else if ((brackets[group_i-1] == ')') && (brackets[group_j+1] == ')'))
      {
         return getPartEnergy(pos_BP_nach, pos_BP_vor, int_seq);
      }

      // free base if left in a BL or IL or ML, pos_BP_vor = closingBP
      else if ((brackets[group_i-1] == '(') && (brackets[group_j+1] == '('))
      {
         //if (BP_Order[pos_BP_nach][3] == 0) ==> BL or IL (left);
         //else ==> ML, closingBP = pos_BP_nach
         // getPartEnergy can be used, since pos_BP_vor and pos_BP_nach are neighbors

         // HERE: pos_BP_nach is predecessor of pos_BP_vor
         return getPartEnergy(pos_BP_vor, pos_BP_nach, int_seq);
      }
      else if (BP_Successors[pos_BP_vor] == BP_Successors[pos_BP_nach])
      {
         //ML
         if (BP_Successors[pos_BP_vor] != numBP)
            return MLEnergy(BP_Successors[pos_BP_vor],int_seq);
         //EL
         else
            return externEnergy(int_seq);
      }

      //forgotten cases
      else
      {
         out("Pos: %d => Forgotten cases!\n", pos_i);
         return 0;
      }
   } // else
}



/*---------------------------------------------------------------------------*/
/* identifies the energy difference that arises by the mutation of a free
   base or a base pair (refers to the whole structure) */
/*---------------------------------------------------------------------------*/

function EnergyDiff(pos_j, sequence)
{
   var pos_i;               // binding pos. of pos. pos_j, if a BP is at pos. pos_j
   var bp_i_new, bp_j_new;  // mutated assignments of the bases of the BP
   var bp_assign_new, base_assign_new; //mutated assignment of the BP and the free base
   var e_old, e_new, e_diff;
   var int_seq, int_seq_new;         //mutated int_seq

   int_seq = char2int(sequence);
   int_seq_new = char2int(sequence);

   //free base
   if (brackets[pos_j] == '.')
   {
      //depending on the location of the base
      for (base_assign_new = 0; base_assign_new < 4; base_assign_new++)
      {
         int_seq_new[pos_j] = base_assign_new;
         e_old = get_BasePart_Energy(pos_j, int_seq);
         e_new = get_BasePart_Energy(pos_j, int_seq_new);
         e_diff = Sub_MinInt(e_old, e_new); // the higher, the better
         /*printf("pos: %d  --  assign_new: %d\n", pos_j, base_assign_new);
         printf("e_old: %f\n", e_old);
         printf("e_new: %f\n", e_new);
         printf("e_dif: %f\n\n", e_diff);*/
         Ediff[pos_j][base_assign_new] = e_diff;
      }
   }
   //BP (only considered at the closing positions)
   else if (brackets[pos_j] == ')')
   {
      pos_i = BP_Order[BP_Pos_Nr[pos_j]][0];
      for (bp_assign_new = 0; bp_assign_new< 6; bp_assign_new++)
      {
         BP2_2(bp_assign_new, bp_i_new, bp_j_new);
         int_seq_new[pos_i] = bp_i_new;
         int_seq_new[pos_j] = bp_j_new;
         e_old = get_BP_Energy(pos_i, int_seq);
         e_new = get_BP_Energy(pos_i, int_seq_new);
         e_diff = Sub_MinInt(e_old, e_new); // the higher, the better
         /*printf("pos_i: %d  --  bp_i_new: %d\n", pos_i, bp_i_new);
         printf("e_old: %f\n", e_old);
         printf("e_new: %f\n", e_new);
         printf("e_dif: %f\n\n", e_diff);*/

         Ediff[pos_j][bp_assign_new] = e_diff;
      }//for (bp_assign)
   }//else if '('
}

/*-------------------------------------------------------------------------*/

                      /* THE LOCAL SEARCH */

           /* some parts used from the Vienna RNA Package */

/*-------------------------------------------------------------------------*/

// function local_search ????


/*-------------------------------------------------------------------------*/

/* shuffle produces a random list by doing len exchanges */
function shuffle(list, length) {
   var i, rnd, temp;
   for (i = 0; i < length; i += 1) {
      rnd = i + Math.floor((Math.random() * (length - i)));
      temp = list[i];
      list[i] = list[rnd];
      list[rnd] = temp;
   }
}

/*-------------------------------------------------------------------------*/

function make_ptable(structure, table) {
   var hx = 0, stack = new Array(structure.length + 1), i, j;
   for (i = 0; i < structure.length; i += 1) {
      switch (structure[i]) {
         case '.' :
            table[i] = -1;
            break;
         case '(' :
            stack[hx++] = i;
            break;
         case ')' :
            j = stack[--hx];
            if (hx < 0) {
               outln("Unbalanced brackets in make_ptable");
            }
            table[i] = j;
            table[j] = i;
            break;
      }
   }
   
   if (hx !== 0) {
      outln("Unbalanced brackets in make_ptable, at end");
   }
}
 
/*-------------------------------------------------------------------------*/

function strncpy(destination, source, num) 
{
   for(var i = 0; i < num; i++)
      destination[i] = source[i];
}

function WALK(i, j) 
{
   strncpy(wstruct, brackets + i, j-i+1);
   strncpy(wstring, string + i, j-i+1);
   if(time_out == 0)
      dist = local_search(wstring, wstruct, i, j, string);
   strncpy(string + i, wstring, j-i+1);
   if((dist > 0) && (give_up)) return -1;
}
 
function inverse_fold(start)
{
   var i, j, jj, o;
   var pt;
   var string, wstring, wstruct, aux;
   var dist=0;
   var precs; //help for identifying the predecessors and successors
   var int_seq;

   start_time = new Date().getTime() / 1000;

   //nc2 = 0;
   j = o = fold_type = 0;

   if (start.length !=struct_len) {
      throw("inverse_fold: start and structure have unequal length");
   }

   string = new Array(struct_len + 1);
   wstring = new Array(struct_len + 1);
   wstruct = new Array(struct_len + 1);
   pt = new Array(struct_len + 1);
   pt[struct_len] = struct_len+1;

   aux = aux_struct(brackets);
   for(i = 0; i < start.length; i++)
      string[i] = start[i];

   make_ptable(brackets, pt);

   precs = GetPrecursors();
   GetSuccessors(precs);
   alloc_Ediff();
   init_Ediff();
   int_seq = char2int(start);

   while (j<struct_len) {
      while ((j<struct_len)&&(brackets[j]!=')')) {
    if (aux[j]=='[') o++;
    if (aux[j]==']') o--;
    j++;
      }
      i=j;
      while (brackets[--i]!='(');  /* doesn't work for open structure */
      if (aux[i]!='[') { i--; j++;}
      while (pt[j]==i) {
    backtrack_type='C';
    if (aux[i]!='[') {
       while (aux[--i]!='[');
       while (aux[++j]!=']');
       /* WALK(i,j); */
    }

    if(WALK(i,j) == -1) {
      backtrack_type = 'F';
      for(i = 0; i < start.length; i++)
         start[i] = string[i];
      zw_time = new Date().getTime() / 1000;
      return dist;
    }
    o--;
    jj = j; i--;
    while (aux[++j]=='.');
    while ((i>=0)&&(aux[i]=='.')) i--;
    if (pt[j]!=i) {
       backtrack_type = (o==0)? 'F' : 'M';
       if (j-jj>8) { 
         if(WALK((i+1),(jj)) == -1) {
            backtrack_type = 'F';
            for(i = 0; i < start.length; i++)
               start[i] = string[i];
            zw_time = new Date().getTime() / 1000;
            return dist;            
         }
       }
       if(WALK((i+1), (j-1)) == -1) {
         backtrack_type = 'F';
         for(i = 0; i < start.length; i++)
            start[i] = string[i];
         zw_time = new Date().getTime() / 1000;
         return dist;
       }
       while ((i>=0) &&(aux[i]==']')) {
          i=pt[i]-1;
          while ((i>=0)&&(aux[i]=='.')) i--;
          if( WALK((i+1), (j-1)) == -1) {
            backtrack_type = 'F';
            for(i = 0; i < start.length; i++)
               start[i] = string[i];
            zw_time = new Date().getTime() / 1000;
            return dist;
          }
       }
    }
      }
   }
}

/*-------------------------------------------------------------------------*/

function inverse_pf_fold(start)
{
   var dist;
   var dang;
   var precs; //help for identifying the predecessors and successors

   start_time = new Date().getTime() / 1000;

   dang=dangles;
   //in the Vienna package dangles is set to 2, but then the energies (evaluated with -d1)
   //don't fit to the probabilies (-d2) (z.B. kommt eine viel hoehere Wsk. fuer eine
   //suboptimale Struktur raus als fuer die mfe-structure)
   //if (dangles!=0) dangles=2;
   if (dangles!=0) dangles=1;
   update_fold_params();    /* make sure there is a valid pair matrix */

   precs = GetPrecursors();
   GetSuccessors(precs);
   alloc_Ediff();
   init_Ediff();

   fold_type=1;
   do_backtrack = 0;

   dist = local_search(start, brackets, 0, struct_len, start);

   dangles=dang;
   return (dist+final_cost);
}

/*---------------------------------------------------------------------------*/

function mfe_cost(string, structure, target)
{
   var energy, distance;

   if (string.length != target.length) {
      throw("unequal length in mfe_cost");
   }

   energy = fold(string, structure);
   distance = bp_distance(target, structure);

   cost2 = energy_of_struct(string, target) - energy;
   return distance;
}

/*---------------------------------------------------------------------------*/
/*****************************************************************
*      return value is no probability but an energy difference   *
*****************************************************************/

function pf_cost(string, structure, target)
{
   var  f, e;

   f = pf_fold(string, structure);
   e = energy_of_struct(string, target);
   return (e-f-final_cost);
}

/*---------------------------------------------------------------------------*/