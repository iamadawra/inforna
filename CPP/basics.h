/***************************************************************************
                       basics.h  -  global variables
                             -------------------
    begin                : 30-08-2004
    copyright            : (C) 2004 by Anke Busch
    email                : anke.busch@informatik.uni-freiburg.de
 ***************************************************************************/

#ifndef _BASICS__
#define _BASICS__

#include <ctype.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include "energy.h"
#include "hairpin_energy.h"
#include "interior_energy.h"
#include "stacking_energy.h"
#include "bulge_energy.h"
#include "multi_energy.h"
#include "end_energy.h"

extern "C"{
             #include <fold.h>
             #include <part_func.h>
             #include <fold_vars.h>
             #include <utils.h>
             #include <inverse.h>
          }

using namespace std;
/**********************************************************************************
*                    Structure- and Type definitions                              *
**********************************************************************************/
typedef unsigned int unt;


/**********************************************************************************
*                                 Makros                                          *
**********************************************************************************/
#define NUM_ACIDS 20
#define NUM_PAM 21


/**********************************************************************************
*                               Constants                                         *
**********************************************************************************/

extern const int MAX_LOOPS;

extern const double T;
extern const double RT;

extern const double MIN_DOUBLE;
extern const double MAX_DOUBLE;
extern const double DOUBLE_DIFF;

extern const int MIN_INT;
extern const int MAX_INT;

/**********************************************************************************
*                         global variabels (extern)                              *
**********************************************************************************/
extern int step;  //is 1 if we are in the initialization step, 2 during the local search (the variable is necessary since constraint mismatches
                  //are allowed only during the second step, but the functions with allow them (PairPenaly, BasePenalty) are used during the
                  //initialization as well

extern char* brackets;
extern int struct_len;
extern char *ElementStruct; 
extern int** BP_Order;
extern int numBP;
extern int* BP_Pos_Nr;
extern double** D;    // Dynamic Recursion Matrix
extern int**** Trace; // Dynamic Traceback Matrix (has 2 dimensions more than D, since in each D-corresponding field the coordinates of all
                      // predecessor fields have to be stored (two coord. per predecessor) and since there could be more than one predecessors
                      // 2 more dimensions are needed)
extern char* iupac_const;
extern int** seq_constraints;
extern int max_mis; //max. number of mismatches that are allowed among the constrained position or in an interval
extern int num_mis; //counter for occuring mismatches
extern int* mis_vec; //vector where for each position is stored whether a mismatch is allowed or not

extern char* best_char_seq; // designed sequence

extern bool free_bases_set2A; // is true, if free bases have to be set to 'A'
extern bool free_bases_set2C;
extern bool free_bases_set2G;
extern bool free_bases_set2U;
extern bool random_init;

extern int search_strategy;       // gives the search strategy ( 1 = adaptive walk, 2 = full local search, 3 = stochastic local search)
extern int neighbour_choice;      // gives the kind of ranking of the neighbors (1 = random, 2 = energy dependent)
extern int only_mutation_is_step; // gives the information how to count the step done during the stochastic local search
extern int step_multiplier;       // maximal number of steps during SLS = allowed_steps * length
extern double p_accept;           // probability to accept worse neighbors during SLS


/**********************************************************************************
*                         Funktionen                                              *
**********************************************************************************/

char int2char(int base);
char* int2char(int* num_seq, int size);
int char2int_base(char base);
int* char2int(char* seq);
int BP2int(int b1, int b2);
void BP2_2(int bp, int &base_i, int &base_j);

double Double_Round(double x);
int Minimum(int a, int b);
double Minimum(double a, double b);
int Maximum(int a, int b);

double* MiniVec(double* vec, int size);
double* MaxiVec(double* vec, int size);
int SumVec(int* vec, int size);

int RandomBase(int number);
int RandomBasePair();
int RandomBasePair(int number);

int Valid_IUPAC(char iu);
int Compare_IUPAC_Base(char iu, int base);
#endif   // _BASICS_
