/**********************************************************
 translates an integerBase to a characterBase (nucleotide)
**********************************************************/

function int2char(base)
{
   if (base == 0) return 'A';
   else if (base == 1) return 'C';
   else if (base == 2) return 'G';
   else if (base == 3) return 'U';
   else
   {
      throw("There's something wrong in your int-base!");
   }
}

/**********************************************************
 translates an integerSeq. to a characterSeq. (nucleotide)
**********************************************************/

function int2char(num_seq, size)
{
   var seq = new Array(size);
   for (var i=0; i<size; i++)
   {
      if (num_seq[i] == 0) seq[i]='A';
      else if (num_seq[i] == 1) seq[i]='C';
      else if (num_seq[i] == 2) seq[i]='G';
      else if (num_seq[i] == 3) seq[i]= 'U';
      else
      {
         throw("There's something wrong in your int-sequence!");
      }
   }
   return seq;
}


/**********************************************************
 translates a characterBase to an integerBase (nucleotide)
**********************************************************/

function char2int_base(base)
{
   if (toupper(base) == 'A') return 0;
   else if (toupper(base) == 'C') return 1;
   else if (toupper(base) == 'G') return 2;
   else if (toupper(base) == 'U') return 3;
   else
   {
      throw("There's something wrong in your char-base!\n");
   }
}

/**********************************************************
 translates a characterSeq. to an integerSeq. (nucleotide)
**********************************************************/

function char2int(seq)
{
   var num_seq = new Array(seq.length);
   for (var i=0; i<seq.length; i++)
      num_seq[i] = char2int_base(seq[i]);
   return num_seq;
}


/**********************************************************
 translates the two bases (integer) of a BP to an integer
 value for the pair (0,1,2,3,4,5)
**********************************************************/

function BP2int(b1, b2)
{
   var kind_of_BP = 0;
   if ( b1==0 && b2==3 ) kind_of_BP = 0;
   else if ( b1==1 && b2==2 ) kind_of_BP = 1;
   else if ( b1==2 && b2==1 ) kind_of_BP = 2;
   else if ( b1==3 && b2==0 ) kind_of_BP = 3;
   else if ( b1==2 && b2==3 ) kind_of_BP = 4;
   else if ( b1==3 && b2==2 ) kind_of_BP = 5;
   else
   {
     throw("Wrong Basepair (" + int2char(b1) + " - " + int2char(b2) + ") !");
   }
   return kind_of_BP;
}


/**********************************************************
 translates a BP (given as an integer) to two integerBases
**********************************************************/

function BP2_2(bp, base_i, base_j)
{
   if (bp == 0) {base_i=0; base_j=3;}
   else if (bp == 1) {base_i=1; base_j=2;}
   else if (bp == 2) {base_i=2; base_j=1;}
   else if (bp == 3) {base_i=3; base_j=0;}
   else if (bp == 4) {base_i=2; base_j=3;}
   else if (bp == 5) {base_i=3; base_j=2;}
   else
   {
      throw("Wrong BasepairNr.!");
   }
}


/**********************************************************
since addition and substraction of two float or double 
values is often not exact enough, calculations using float
or double values were rounded to 3 positions after decimal
point: the number is multiplied by 1000, 0.5 is added (if
the number is negative, 0.5 is substracted), then turn the
number to an integer, i.e. the remaining positions after
decimal point are cut. 
(The procedure is similar to rounding but 0.5 is added
before.)
**********************************************************/

function Double_Round(x)
{
   var x_int;
   if (x>=0)
      x_int = Math.round((x*1000)+0.5);
   else
      x_int = Math.round((x*1000)-0.5);

   return x_int/1000;
}



/************************************************************
 finds the minimum of two integers
************************************************************/

function Minimum(a, b)
{
   if (a<b)
      return a;
   else
      return b;
}


/************************************************************
 finds the maximum of two integers
************************************************************/

function Maximum(a, b)
{
   if (a>b)
      return a;
   else
      return b;
}


/************************************************************
 sums up all values of a vector
************************************************************/

function SumVec(vec, size)
{
   var sum = 0;
   for (var i=0; i<size; i++)
      sum += vec[i];
   return sum;
}

/************************************************************
 finds the minimal value in a vector (returns the value and
 its coordinate)
************************************************************/

function MiniVec(vec, size)
{
   var min = new Array(2);
   min[0] = 1000000.0; //coord.
   min[1] = 1000000.0; //value
   for ( var i=0; i<size; i++)
      if (vec[i] < min[1])
      {
         min[1] = vec[i];
         min[0] = i;
      }
   return min;
}

/************************************************************
 finds the maximal value in a vector (returns the value and
 its coordinate)
************************************************************/

function MaxiVec(vec, size)
{
   var max = new Array(2);
   max[0] = -1000000; //coord.
   max[1] = -1000000; //value
   for (var i=0; i<size; i++)
      if (vec[i] > max[1])
      {
         max[1] = vec[i];
         max[0] = i;
      }
   return max;
}


/*************************************************************
 gives a random base (depending on the size of the set of 
 valid bases)
*************************************************************/

function RandomBase(int number)
{
   var zufall;
   zufall = Math.random();

   if (number == 4)
   {
      if (zufall <= 0.25)
         return 0;
      else if ((zufall > 0.25) && (zufall <= 0.5))
         return 1;
      else if ((zufall > 0.5) && (zufall <= 0.75))
         return 2;
      else
         return 3;
   }
   else if (number == 3)
   {
      if (zufall <= 0.3333)
         return 0;
      else if ((zufall > 0.3333) && (zufall <= 0.6667))
         return 1;
      else
         return 2;
   }
   else if (number == 2)
   {
      if (zufall <= 0.5)
         return 0;
      else
         return 1;
   }
   else
      return 0;
}


/*************************************************************
 finds a random base pair
*************************************************************/
function RandomBasePair()
{
   var zufall;
   zufall = Math.random();

   if (zufall <= 0.1667)
      return 0;
   else if ((zufall > 0.1667) && (zufall <= 0.3333))
      return 1;
   else if ((zufall > 0.3333) && (zufall <= 0.5))
      return 2;
   else if ((zufall > 0.5) && (zufall <= 0.6667))
      return 3;
   else if ((zufall > 0.6667) && (zufall <= 0.8333))
      return 4;
   else 
      return 5;
}


/*************************************************************
 finds a random base pair depending on constraints
*************************************************************/

function RandomBasePair(number)
{
   var zufall;
   zufall = Math.random();

   if (number == 6)
   {
      if (zufall <= 0.1667)
         return 0;
      else if ((zufall > 0.1667) && (zufall <= 0.3333))
         return 1;
      else if ((zufall > 0.3333) && (zufall <= 0.5))
         return 2;
      else if ((zufall > 0.5) && (zufall <= 0.6667))
         return 3;
      else if ((zufall > 0.6667) && (zufall <= 0.8333))
         return 4;
      else 
         return 5;
   }
   else if (number == 5)
   {
      if (zufall <= 0.2)
         return 0;
      else if ((zufall > 0.2) && (zufall <= 0.4))
         return 1;
      else if ((zufall > 0.4) && (zufall <= 0.6))
         return 2;
      else if ((zufall > 0.6) && (zufall <= 0.8))
         return 3;
      else
         return 4;
   }
   else if (number == 4)
   {
      if (zufall <= 0.25)
         return 0;
      else if ((zufall > 0.25) && (zufall <= 0.5))
         return 1;
      else if ((zufall > 0.5) && (zufall <= 0.75))
         return 2;
      else
         return 3;
   }
   else if (number == 3)
   {
      if (zufall <= 0.3333)
         return 0;
      else if ((zufall > 0.3333) && (zufall <= 0.6667))
         return 1;
      else
         return 2;
   }
   else if (number == 2)
   {
      if (zufall <= 0.5)
         return 0;
      else
         return 1;
   }
   else
      return 0;
}


/*************************************************************
 test, whether a base is valid in IUPAC-code
*************************************************************/
function Valid_IUPAC(iu)
{
   if ((iu == 'A') || (iu == 'C') || (iu == 'G') || (iu == 'U') || (iu == 'M') || (iu == 'R') || (iu == 'W') || (iu == 'S') || (iu == 'Y') || (iu == 'K') || (iu == 'V') || (iu == 'H') || (iu == 'D') || (iu == 'B') || (iu == 'N'))
      return 1;
   else
      return 0;
}

/*************************************************************
 test, whethter a base is valid for a given UIPAC-Code base
*************************************************************/
function Compare_IUPAC_Base(iu, base)
{
   if (base == 0)
   {
      if ((iu == 'A') || (iu == 'M') || (iu == 'R') || (iu == 'W') || (iu == 'V') || (iu == 'H') || (iu == 'D') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else if (base == 1)
   {
      if ((iu == 'C') || (iu == 'M') || (iu == 'S') || (iu == 'Y') || (iu == 'V') || (iu == 'H') || (iu == 'B') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else if (base == 2)
   {
      if ((iu == 'G') || (iu == 'R') || (iu == 'S') || (iu == 'K') || (iu == 'V') || (iu == 'D') || (iu == 'B') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else if (base == 3)
   {
      if ((iu == 'U') || (iu == 'W') || (iu == 'Y') || (iu == 'K') || (iu == 'H') || (iu == 'D') || (iu == 'B') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else
   {
      throw("Wrong Base!");
   }
}





