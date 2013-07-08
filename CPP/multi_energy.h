#ifndef _MULTI_ENERGY__
#define _MULTI_ENERGY__

#include <stdlib.h>
#include "basics.h"
#include "constraints.h"

using namespace std;

int* MultiLoopConnections(int stem_num, vector<int> & stem_ends, int order_pos_of_closing_bp, vector< vector<int> > & BaseConnections, vector< vector<int> > & BasePairConnections);
void MultiLoopConnections(int* stem_ends, int order_pos_of_closing_bp, int** BaseConnections, int *BaseSizes, int** BasePairConnections, int* BasePairSizes);

double BestConnectionEnergy(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection);
double ConnectionEnergy(const int* base_connection, int base_size, const int* pair_connection, int bp_size, int* int_seq);

void MLBestEnergy(int bp_pos, vector<int> & stem_ends);
double MLEnergy(int bp_pos, int* int_seq);
int* FindStemEnds(int closingBP);

int* ConnectionBestFreeBases(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection);

#endif   // _MULTI_ENERGY_
