#ifndef _SEARCH__
#define _SEARCH__

#include <stdlib.h>
#include <limits.h>
#include "basics.h"
#include "constraints.h"
#include "struct.h"


using namespace std;


float inverse_fold(char *start);
/* find sequences with predifined structure.
   the found sequence is written to start,
   return value is
      energy_of_struct(start, target) - fold(start, structure),
   i.e. 0. if search was successful; */

float inverse_pf_fold(char *start);
/*  inverse folding maximising the frequency of target in the
    ensemble of structures, final sequence is written to start, returns
       energy_of_struct(start, target) - part_func(start, structure)
*/


void Pos2BP_Pos();
int** GetPrecursors();
void GetSuccessors(int** precs);
void alloc_Ediff();
void init_Ediff();
void print_Ediff();
void print_av_Ediff();
void make_av_Ediff();
void print_max_Ediff();
void make_max_Ediff();

int compare_Pos(const void *a, const void *b);
int* make_mut_pos_list();

double getPartEnergy(int bp_pos, int bp_pos_before, int* int_seq);
double get_BP_Energy(int pos_i, int* int_seq);
double get_BasePart_Energy(int pos_i, int* int_seq);

void EnergyDiff(int pos_i, char* sequence);

double local_search(char *start, char *target, int pos_i, int pos_j, char* whole_seq);
void   shuffle(int *list, int len);
void   make_ptable(char *structure, int *table);
double  mfe_cost(char *, char*, char *);
double  pf_cost(char *, char *, char *);



#endif   // _SEARCH_

