
#include "struct.h"

/*--------------------------------------------------------------------------*/
/************************************************
checks, whether the user given structure is valid

the following structure should be valid as well,
since there are no parameters that qualify one
sequence more than another: "..........."

copied from: Vienna Package + extended:
************************************************/

int check_brackets(char *line)
{
  int i,o,bonds;

  i=o=bonds=0;
  while( line[i] ){
    switch(line[i]) {
    case '(' :
      o++;
      i++;
      bonds++;
      break;
    case '.' :
      i++;
      break;
    case ')' : 
      i++;
      o--;
      if(o<0) return 0;
      break;
    default:
      return 0;
    }
  }
  if (o>0) return 0;
  if (bonds == 0) return 0;
  return 1;
}

/*--------------------------------------------------------------------------*/

/*****************************************************
for each base the pairing base is stored,
for free bases -1 is stored

copied from: Vienna Package
*****************************************************/

int* make_BasePair_Table(char *structure)
{
   int j,hx;
   int *stack, *table;

   hx=0;
   stack = (int *) malloc(sizeof(int)*(strlen(structure)+1));
   table = (int *) malloc(sizeof(int)*(strlen(structure)+1));

   for (unsigned int i=0; i<strlen(structure); i++)
   {
      switch (structure[i])
      {
          case '.':
                     table[i]= -1;
                     break;
          case '(':
                     stack[hx++]=i;
                     break;
          case ')':
                     j = stack[--hx];
                     if (hx<0)
                     {
                        cerr << structure << endl;
                        cerr << "unbalanced brackets in make_BasePair_Table";
                     }
                     table[i]=j;
                     table[j]=i;
                     break;
      }
   }
   if (hx!=0)
   {
      cerr << structure << endl;
      cerr << "unbalanced brackets in make_BasePair_Table";
   }
   free(stack);
   return table;
}

/*-------------------------------------------------------------------------*/
/*****************************************
counts the number of BPs
*****************************************/

int NumOfBP(char* structure)
{
   int num = 0;
   for (unsigned int i=0; i<strlen(structure); i++)
      if (structure[i] == '(')
         num++;
   return num;
}

/*-------------------------------------------------------------------------*/
/**************************************************
translates the structure such that the last BP of
a stem is represented by []

copied from: Vienna Package
***************************************************/

char *aux_struct(const char* structure )
{  
   short        *match_paren;
   int          i, o, p;
   char         *string;
   
   string = (char *) malloc(sizeof(char)*(strlen(structure)+1));
   match_paren = (short *) malloc(sizeof(short)*(strlen(structure)/2+1));
   strcpy(string, structure);

   i = o = 0;
   while (string[i])
   {
      switch (string[i])
      {
          case '.':  break;
          case '(':
                     match_paren[++o]=i;
                     break;
          case ')':
                     p=i;
                     while ((string[p+1]==')')&&(match_paren[o-1]==match_paren[o]-1))
                     {
                        p++;
                        o--;
                     }
                     string[p]=']';
                     i=p;
                     string[match_paren[o]]='[';
                     o--;
                     break;
          default:
                     cerr << "Junk in structure at aux_structure\n";
      }
      i++;
   }
   free(match_paren);
   return(string);
}

/*---------------------------------------------------------------------------*/

/**********************************************************
translates the usual bracket notation of the structure (or the
structure after aux_struct) in a representation that just gives
structual elements and their sizes,
e.g.: ((((((H3)S3)((H3)S3)M4)S2)E2)

copied from: Vienna Package
**********************************************************/

char *right_Element_Structure(const char *structure )
{
   short *bulge, *loop;

   int            i, lp, p, l;
   unsigned int   k;
   char           *string, *Shapiro, *temp, tt[10];

   int    loop_size[MAX_LOOPS];       /* contains loop sizes of a structure */
   int    helix_size[MAX_LOOPS];      /* contains helix sizes of a structure */
   int    loop_degree[MAX_LOOPS];     /* contains loop degrees of a structure */
   int    loops;                  /* n of loops and stacks in a structure */
   int    unpaired, pairs;        /* n of unpaired digits and pairs */

   bulge = (short *) malloc(sizeof(short)*(strlen(structure)/3+1));
   loop = (short *) malloc(sizeof(short)*(strlen(structure)/3+1));
   temp = (char *) malloc(4*strlen(structure)+1);

   for (i = 0; i < MAX_LOOPS; i++)
      loop_size[i] = helix_size[i] = 0;

   loop_degree[0]=0;         /* open structure has degree 0 */
   pairs = unpaired = loops = lp = 0;
   loop[0]=0;

   string = aux_struct( structure );

   i=p=l=0;
   while (string[i])
   {
      switch(string[i])
      {
         case '.':
                    unpaired++;
                    loop_size[loop[lp]]++;
                    break;
         case '[':
                    temp[l++]='(';
                    temp[l++]='(';
                    if ((i>0)&&(string[i-1]=='('))
                       bulge[lp]=1;
                    lp++;
                    loop_degree[++loops]=1;
                    loop[lp]=loops;
                    bulge[lp]=0;
                    break;
         case ')':
                    if (string[i-1]==']')
                       bulge[lp]=1;
                    p++;
                    break;
         case ']':
                    if (string[i-1]==']')
                       bulge[lp]=1;
                    switch (loop_degree[loop[lp]])
                    {
                       case 1:  temp[l++]='H';
                                break;           /* hairpin */
                       case 2:
                                if (bulge[lp]==1)
                                   temp[l++] = 'B';             /* bulge */
                                else
                                   temp[l++] = 'I';             /* internal loop */
                                break;
                       default: temp[l++] = 'M';                /* multiloop */
                    }
                    helix_size[loop[lp]]=p+1;

                    sprintf(tt, "%d)" , loop_size[loop[lp]]);
                    for(k=0; k<strlen(tt); k++)
                       temp[l++] = tt[k];
                    sprintf(tt, "S%d)" , helix_size[loop[lp]]);
                    for(k=0; k<strlen(tt); k++)
                       temp[l++] = tt[k];

                    pairs+=p+1;
                    p=0;
                    loop_degree[loop[--lp]]++;
                    break;
      }
      i++;
   }

   *tt = '\0';
   if (loop_size[0])
      sprintf(tt, "E%d)" , loop_size[0]);

   for(k=0; k<strlen(tt); k++)
      temp[l++] = tt[k];
   temp[l]='\0';

   if (loop_size[0])
      l++;

   Shapiro = (char *) malloc(sizeof(char)*(l+1));
   if (loop_size[0])
   {
      Shapiro[0]='(';
      strcpy(Shapiro+1, temp);
   }
   else
      strcpy(Shapiro, temp);

   free(string);
   free(temp);
   free(loop); free(bulge);

   return Shapiro;
}




/*---------------------------------------------------------------------------*/


/**********************************************************

translates the structure representation in another that
shows the element and its size at both: the opening and the
closing bracket

***********************************************************/

char* Element_Structure(char* right_structure)
{
   int k;
   char* newStruct;
   vector<char> tmp;              //helping vector, contains all closing brackets and its elements
   vector<char> newStructInverse; //the inverse of the new structure, since the structure is checked from tail to head

   for (int i=strlen(right_structure)-1; i>=0; i--)
   {
      if (right_structure[i] != '(')
      {
         tmp.push_back(right_structure[i]);
         newStructInverse.push_back(right_structure[i]);
      }
      else
      {
         k=0;
         while (tmp[tmp.size()-1-k] != ')')
            k++;
         // k is the number of entries of tmp that are written to newStructInv.
         for (int j=k-1; j>=0; j--)
            newStructInverse.push_back(tmp[tmp.size()-1-j]);
         newStructInverse.push_back('(');
         for (int j=0; j<=k; j++)
            tmp.pop_back();
      }
   }

   newStruct = (char*) malloc(sizeof(char)*newStructInverse.size());
   for (unsigned int i=0; i<newStructInverse.size(); i++)
      newStruct[i] = newStructInverse[newStructInverse.size()-i-1];

   return newStruct;
}


/*---------------------------------------------------------------------------*/


/**********************************************************
if there are more than one structure fragments connected by
some free bases without a closing BP (= external loop)

This external loop may contain dangling ends and its energy
values are stored in the last row of D and Trace. Thus
an equivalent position in BP_Order is needed and there,
the number of free bases and stems is stored.
***********************************************************/

int* ClosingStructure(int* bpTable)
{
   int* stems_and_freeBases;
   stems_and_freeBases = (int*) malloc(sizeof(int)*2);
   int stems = 0;
   int freeB = 0;
   int i = 0;

   while (i < (int)strlen(brackets))
   {
      if (brackets[i] == '.')
      {
         freeB++;
         i++;
      }
      else if (brackets[i] == '(')
      {
         stems++;
         i = bpTable[i] + 1;
      }
   }
   stems_and_freeBases[0] = freeB;
   stems_and_freeBases[1] = stems;
   return stems_and_freeBases;
}


/*---------------------------------------------------------------------------*/




