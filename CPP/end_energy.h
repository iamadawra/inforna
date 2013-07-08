#ifndef _END_ENERGY__
#define _END_ENERGY__

#include <stdlib.h>
#include "basics.h"
#include "struct.h"
#include "constraints.h"

using namespace std;

void EndConnections(int stem_num, int* vorgaenger, vector< vector<int> > & EndBaseConnections, vector< vector<int> > & EndBasePairConnections);
void EndConnections(int* vorgaenger, int** EndBaseConnections, int *EndBaseSizes, int** EndBasePairConnections, int *EndBasePairSizes);

double BestEndConnectionEnergy(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection);
double EndConnectionEnergy(const int* base_connection, int base_size, const int* pair_connection, int pair_size, int* int_seq);

void externBestEnergy();
double externEnergy(int* int_seq);

int* EndConnectionBestFreeBases(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection);

#endif   // _END_ENERGY_
