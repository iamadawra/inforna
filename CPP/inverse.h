#ifndef _INVERSE__
#define _INVERSE__

#include <stdlib.h>
#include <limits.h>
#include "basics.h"
#include "struct.h"
#include "constraints.h"

using namespace std;

void getOrder(char* structure);
int* StemsInML(int first_pos, int BP_start_pos);
float Recursion();
int* Traceback();
double Random_Init();

#endif   // _INVERSE_

