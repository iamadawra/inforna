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