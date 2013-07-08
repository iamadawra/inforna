#ifndef _STRUCT__
#define _STRUCT__

#include <stdlib.h>
#include "basics.h"
#include "constraints.h"

using namespace std;

int check_brackets(char *line);
int* make_BasePair_Table(char *structure);
int NumOfBP(char* structure);
char *aux_struct(const char* structure );
char *right_Element_Structure(const char *structure );
char* Element_Structure(char* right_structure);
int* ClosingStructure(int* bpTable);


#endif   // _STRUCT_
