var MAXALPHA = 20;
var TIME_OUT_TIME = 3600;

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