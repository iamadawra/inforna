#include "basics.h"

/**********************************************************
 translates an integerBase to a characterBase (nucleotide)
**********************************************************/

char int2char(int base)
{
   if (base == 0) return 'A';
   else if (base == 1) return 'C';
   else if (base == 2) return 'G';
   else if (base == 3) return 'U';
   else
   {
      cerr << "There's something wrong in your int-base!";
      exit(1);
   }
}

/**********************************************************
 translates an integerSeq. to a characterSeq. (nucleotide)
**********************************************************/

char* int2char(int* num_seq, int size)
{
   char* seq = (char*) malloc(sizeof(char)*size);
   for (int i=0; i<size; i++)
   {
      if (num_seq[i] == 0) seq[i]='A';
      else if (num_seq[i] == 1) seq[i]='C';
      else if (num_seq[i] == 2) seq[i]='G';
      else if (num_seq[i] == 3) seq[i]= 'U';
      else
      {
         cerr << "There's something wrong in your int-sequence!";
         exit(1);
      }
   }
   return seq;
}


/**********************************************************
 translates a characterBase to an integerBase (nucleotide)
**********************************************************/

int char2int_base(char base)
{
   if (toupper(base) == 'A') return 0;
   else if (toupper(base) == 'C') return 1;
   else if (toupper(base) == 'G') return 2;
   else if (toupper(base) == 'U') return 3;
   else
   {
      printf("There's something wrong in your char-base!\n");
      exit(1);
   }
}

/**********************************************************
 translates a characterSeq. to an integerSeq. (nucleotide)
**********************************************************/

int* char2int(char* seq)
{
   int* num_seq = (int*) malloc(sizeof(int)*strlen(seq));
   for (int i=0; i<(int)strlen(seq); i++)
      num_seq[i] = char2int_base(seq[i]);
   return num_seq;
}


/**********************************************************
 translates the two bases (integer) of a BP to an integer
 value for the pair (0,1,2,3,4,5)
**********************************************************/

int BP2int(int b1, int b2)
{
   int kind_of_BP = 0;
   if ( b1==0 && b2==3 ) kind_of_BP = 0;
   else if ( b1==1 && b2==2 ) kind_of_BP = 1;
   else if ( b1==2 && b2==1 ) kind_of_BP = 2;
   else if ( b1==3 && b2==0 ) kind_of_BP = 3;
   else if ( b1==2 && b2==3 ) kind_of_BP = 4;
   else if ( b1==3 && b2==2 ) kind_of_BP = 5;
   else
   {
      cerr << "Wrong Basepair (" << int2char(b1) << " - " << int2char(b2) << ") !" << endl;
      exit(1);
   }
   return kind_of_BP;
}


/**********************************************************
 translates a BP (given as an integer) to two integerBases
**********************************************************/

void BP2_2(int bp, int &base_i, int &base_j)
{
   if (bp == 0) {base_i=0; base_j=3;}
   else if (bp == 1) {base_i=1; base_j=2;}
   else if (bp == 2) {base_i=2; base_j=1;}
   else if (bp == 3) {base_i=3; base_j=0;}
   else if (bp == 4) {base_i=2; base_j=3;}
   else if (bp == 5) {base_i=3; base_j=2;}
   else
   {
      cerr << "Wrong BasepairNr.!" << endl;
      exit(1);
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

double Double_Round(double x)
{
   int x_int;
   if (x>=0)
      x_int = (int)((x*1000)+0.5);
   else
      x_int = (int)((x*1000)-0.5);

   return (double) x_int/1000;
}



/************************************************************
 finds the minimum of two integers
************************************************************/

int Minimum(int a, int b)
{
   if (a<b)
      return a;
   else
      return b;
}

/************************************************************
 finds the minimum of two double values
************************************************************/

double Minimum(double a, double b)
{
   if (a<b)
      return a;
   else
      return b;
}


/************************************************************
 finds the maximum of two integers
************************************************************/

int Maximum(int a, int b)
{
   if (a>b)
      return a;
   else
      return b;
}


/************************************************************
 sums up all values of a vector
************************************************************/

int SumVec(int* vec, int size)
{
   int sum = 0;
   for (int i=0; i<size; i++)
      sum += vec[i];
   return sum;
}

/************************************************************
 finds the minimal value in a vector (returns the value and
 its coordinate)
************************************************************/

double* MiniVec(double* vec, int size)
{
   double* min = (double*) malloc(sizeof(double)*2);
   min[0] = MAX_DOUBLE; //coord.
   min[1] = MAX_DOUBLE; //value
   for ( int i=0; i<size; i++)
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

double* MaxiVec(double* vec, int size)
{
   double* max = (double*) malloc(sizeof(double)*2);
   max[0] = MIN_INT; //coord.
   max[1] = MIN_INT; //value
   for (int i=0; i<size; i++)
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

int RandomBase(int number)
{
   double zufall;
   zufall = rand();
   zufall /= RAND_MAX;

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
int RandomBasePair()
{
   double zufall;
   zufall = rand();
   zufall /= RAND_MAX;

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

int RandomBasePair(int number)
{
   double zufall;
   zufall = rand();
   zufall /= RAND_MAX;
   cout << "zufall: " << zufall << endl;

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
int Valid_IUPAC(char iu)
{
   if ((iu == 'A') || (iu == 'C') || (iu == 'G') || (iu == 'U') || (iu == 'M') || (iu == 'R') || (iu == 'W') || (iu == 'S') || (iu == 'Y') || (iu == 'K') || (iu == 'V') || (iu == 'H') || (iu == 'D') || (iu == 'B') || (iu == 'N'))
      return 1;
   else
      return 0;
}

/*************************************************************
 test, whethter a base is valid for a given UIPAC-Code base
*************************************************************/
int Compare_IUPAC_Base(char iu, int base)
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
      cerr << "Wrong Base!" << endl;
      exit(1);
   }
}





