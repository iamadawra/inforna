//Needs Vector Functionality
//STRCOPY IN JS
//VECTOR LIBRARY DEPENDENCIES
//What the hell is \0 doing in js??!

/*--------------------------------------------------------------------------*/
/************************************************
checks, whether the user given structure is valid

the following structure should be valid as well,
since there are no parameters that qualify one
sequence more than another: "..........."

copied from: Vienna Package + extended:
************************************************/

function check_brackets(line)
{
  var i,o,bonds;

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

function make_BasePair_Table(structure)
{
   var j,hx;
   var stack, table;

   hx=0;
   stack = new Array((structure.length)+1);
   table = new Array((structure.length)+1);

   for (var i=0; i<(structure.length); i++)
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
                        throw("unbalanced brackets in make_BasePair_Table");
                     }
                     table[i]=j;
                     table[j]=i;
                     break;
      }
   }
   if (hx!=0)
   {
      throw("unbalanced brackets in make_BasePair_Table");
   }
   return table;
}

/*-------------------------------------------------------------------------*/
/*****************************************
counts the number of BPs
*****************************************/

function NumOfBP(structure)
{
   var num = 0;
   for (var i=0; i<(structure.length); i++)
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

function aux_struct(structure)
{  
   var        match_paren;
   var          i, o, p;
   var         string;
   
   string = new Array((structure.length)+1);
   match_paren = new Array((structure.length)/2+1);
   strcpy(string, structure); //STRCOPY IN JS?????????

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
                     throw("Junk in structure at aux_structure\n");
      }
      i++;
   }
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

function right_Element_Structure(structure)
{
   var bulge, loop;

   var            i, lp, p, l;
   var            k;
   var            string, Shapiro, temp, tt[10];

   MAX_LOOPS = 2000;

   var    loop_size[MAX_LOOPS];       /* contains loop sizes of a structure */
   var    helix_size[MAX_LOOPS];      /* contains helix sizes of a structure */
   var    loop_degree[MAX_LOOPS];     /* contains loop degrees of a structure */
   var    loops;                  /* n of loops and stacks in a structure */
   var    unpaired, pairs;        /* n of unpaired digits and pairs */

   bulge = new Array((structure.length)/3+1);
   loop = new Array((structure.length)/3+1);
   temp = new Array(4*(structure.length)+1);

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
                    for(k=0; k<strlen(tt); k++)
                       temp[l++] = tt[k];
                    for(k=0; k<strlen(tt); k++)
                       temp[l++] = tt[k];

                    pairs+=p+1;
                    p=0;
                    loop_degree[loop[--lp]]++;
                    break;
      }
      i++;
   }

   tt = '\0';

   for(k=0; k<strlen(tt); k++)
      temp[l++] = tt[k];
   temp[l]='\0';

   if (loop_size[0])
      l++;

   Shapiro = new Array(l+1);
   if (loop_size[0])
   {
      Shapiro[0]='(';
      strcpy(Shapiro+1, temp);   //STRCOPY???????????????
   }
   else
      strcpy(Shapiro, temp);     //STRCOPY???????????????

   return Shapiro;
}




/*---------------------------------------------------------------------------*/


/**********************************************************

translates the structure representation in another that
shows the element and its size at both: the opening and the
closing bracket

***********************************************************/

function Element_Structure(right_structure)
{
   var k;
   var newStruct;     // VECTORS ????????????????????????????????????
   vector<char> tmp;              //helping vector, contains all closing brackets and its elements
   vector<char> newStructInverse; //the inverse of the new structure, since the structure is checked from tail to head

   for (var i=strlen(right_structure)-1; i>=0; i--)
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
         for (var j=k-1; j>=0; j--)
            newStructInverse.push_back(tmp[tmp.size()-1-j]);
         newStructInverse.push_back('(');
         for (var j=0; j<=k; j++)
            tmp.pop_back();
      }
   }

   newStruct = new Array(newStructInverse.size());
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

function ClosingStructure(bpTable)
{
   var stems_and_freeBases;
   stems_and_freeBases = new Array(2);
   var stems = 0;
   var freeB = 0;
   var i = 0;

   while (i < (brackets.length))
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




