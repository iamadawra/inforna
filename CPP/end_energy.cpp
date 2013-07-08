#include "end_energy.h"

/*******************************************************
********************************************************

            DEALS WITH THE EXTERNAL LOOP

********************************************************
*******************************************************/

/*******************************************************
 identifies connections concerning the closing base pairs
 of the stems and concerning the free bases adjacent to them

 stem_num = number of stems
 vorgaenger = contains all stem ends and their positions in BP_Order
 BaseConnections = contains row-wise the connections concerning bases
 BasePairConnections = contains row-wise the connections concerning base pairs
*******************************************************/

void EndConnections(int stem_num, int* vorgaenger, vector< vector<int> > & EndBaseConnections, vector< vector<int> > & EndBasePairConnections)
{
   int* bases_between_stems = (int*) malloc(sizeof(int)*stem_num);  //contains the number of bases between the stems (in each coord. the bases "before" the
                                                                    //stem, only at the last coord. of bases_between_stems the free bases after
                                                                    //the last stem and before the closing BP are stored
   vector<int> one_connection;        //help vector for one connection

   // number of free bases between the stem-ending base pairs
   for (int i=0; i<stem_num-1; i++)
      bases_between_stems[i] = BP_Order[vorgaenger[i+1]][0] - BP_Order[vorgaenger[i]][1] - 1;
   bases_between_stems[stem_num-1] = (int)strlen(brackets) - BP_Order[vorgaenger[stem_num-1]][1] - 1;

   // finding all connections: (concerning free bases)
   //**************************************************
   for (int i=0; i<stem_num-1; i++)
   {
      // if there is only one free base between the stems and if there is at least one free base after the next stem
      if (bases_between_stems[i] >= 1)
      {
         one_connection.push_back(BP_Order[vorgaenger[i]][1]+1);
         if (bases_between_stems[i] > 1)
         {
            EndBaseConnections.push_back(one_connection);
            one_connection.clear();
            one_connection.push_back(BP_Order[vorgaenger[i+1]][0]-1);
         }
      }
      else
      {
         EndBaseConnections.push_back(one_connection);
         one_connection.clear();
      }
   }

   if (bases_between_stems[stem_num-1] >= 1)
   {
      one_connection.push_back(BP_Order[vorgaenger[stem_num-1]][1]+1);
      EndBaseConnections.push_back(one_connection);
      one_connection.clear();
   }
   else
   {
      EndBaseConnections.push_back(one_connection);
      one_connection.clear();
   }


   // finding all connections: (concerning stems)
   //*************************************************
   one_connection.push_back(vorgaenger[0]);

   for (int i=0; i<stem_num-1; i++)
   {
      if (bases_between_stems[i] == 1)
         one_connection.push_back(vorgaenger[i+1]);
      else
      {
         EndBasePairConnections.push_back(one_connection);
         one_connection.clear();
         one_connection.push_back(vorgaenger[i+1]);
      }
   }
   EndBasePairConnections.push_back(one_connection);
   one_connection.clear();
}


/*******************************************************
 identifies connections concerning the closing base pairs
 of the stems and concerning the free bases adjacent to them
 (similar to the other EndConnections-function, but without 
 using C++ vectors)

 vorgaenger = contains all stem ends of leaving stems and their positions in BP_Order
 BaseConnections = contains row-wise the connections concerning bases
 BasePairConnections = contains row-wise the connections concerning base pairs
 EndBaseSizes = sizes of the EndBaseConnections
 EndBasePairSizes = sizes of the EndBasePairConnections
*******************************************************/

void EndConnections(int* vorgaenger, int** EndBaseConnections, int *EndBaseSizes, int** EndBasePairConnections, int *EndBasePairSizes)
{
   int stem_num = Maximum(BP_Order[numBP][3],1);
   int num_bp_con = 0;    //current number of BasePairConnections
   int num_base_con = 0;  //current number of BaseConnections
   int bp, base; //current position in the current connection
   int* bases_between_stems = (int*) malloc(sizeof(int)*stem_num);  //contains the number of bases between the stems (in each coord. the bases "before" the
                                                                    //stem, only at the last coord. of bases_between_stems the free bases after
                                                                    //the last stem and before the closing BP are stored

   // number of free bases between the stem-ending base pairs
   for (int i=0; i<stem_num-1; i++)
      bases_between_stems[i] = BP_Order[vorgaenger[i+1]][0] - BP_Order[vorgaenger[i]][1] - 1;
   bases_between_stems[stem_num-1] = struct_len - BP_Order[vorgaenger[stem_num-1]][1] - 1;

   // finding all connections: (concerning free bases)
   //**************************************************
   bp = 0;
   base = 0;

   for (int i=0; i<stem_num-1; i++)
   {
      // if there is only one free base between the stems and if there is at least one free base after the next stem
      if (bases_between_stems[i] >= 1)
      {
         EndBaseConnections[num_base_con][base++] = BP_Order[vorgaenger[i]][1]+1;
         if (bases_between_stems[i] > 1)
         {
            EndBaseSizes[num_base_con] = base;
            num_base_con++;
            base = 0;
            EndBaseConnections[num_base_con][base++] = BP_Order[vorgaenger[i+1]][0]-1;
         }
      }
      else
      {
         EndBaseSizes[num_base_con] = base;
         num_base_con++;
         base = 0;
      }
   }

   if (bases_between_stems[stem_num-1] >= 1)
   {
      EndBaseConnections[num_base_con][base++] = BP_Order[vorgaenger[stem_num-1]][1]+1;
      EndBaseSizes[num_base_con] = base;
      num_base_con++;
      base = 0;
   }
   else
   {
      EndBaseSizes[num_base_con] = base;
      num_base_con++;
      base = 0;
   }


   // finding all connections: (concerning stems)
   //*************************************************
   EndBasePairConnections[num_bp_con][bp++] = vorgaenger[0];

   for (int i=0; i<stem_num-1; i++)
   {
      if (bases_between_stems[i] == 1)
         EndBasePairConnections[num_bp_con][bp++] = vorgaenger[i+1];
      else
      {
         EndBasePairSizes[num_bp_con] = bp;
         num_bp_con++;
         bp = 0;
         EndBasePairConnections[num_bp_con][bp++] = vorgaenger[i+1];
      }
   }
   EndBasePairSizes[num_bp_con] = bp;
   num_bp_con++;
   bp = 0;
}



/******************************************************
gives the MINIMAL (=best) free energy for a connection of the
external loop if the base pair assignments of the stem-ending
base pairs are known

base_connection = a connection concerning free bases
pair_connection = a connection concerning stem-ending base pairs
bp_at_pair_connection = a BP assignment of pair_connection
******************************************************/

double BestEndConnectionEnergy(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection)
{
   double energy = 0.0;
   double min = MAX_DOUBLE;
   int base_size = (int)base_connection.size();
   int pair_size = (int)pair_connection.size();
   double base_energy;

   for (int i=0; i<base_size; i++)
   {
      min = MAX_DOUBLE;

      /*There are 4 cases:  (a) BP b BP b
                            (b) b BP b BP b
                            (c) b BP b BP
                            (d) BP b BP*/

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         // corresponds to cases (b) and (c) (from the second free base on)
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               //energy and penalty
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }
      }

      if (i < pair_size)
      {
         // usual case: base is located left of (= before) the stem
         // corresponds to cases (b) and (c) (free bases in front of a following stem)
         if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located right of (= after) the stem
         // corresponds to cases (a) and (d) (free bases after the stem)
         if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }
      }

      if ( i+1 < pair_size)
      {
         // usual case: base is located left of (=before) the stem
         // correponds to cases (a) and (d) (free base in front of the stem)
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][1]+4*bp_at_pair_connection[i+1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }
      }
      energy = Sum_MaxDouble(energy,min);
   }
   return energy;
}




/******************************************************
gives the free energy for a connection of the
external loop if the base pair assignments of the stem-ending
base pairs and the assignments of the free bases are known

base_connection = a connection concerning free bases
pair_connection = a connection concerning stem-ending base pairs
******************************************************/

double EndConnectionEnergy(const int* base_connection, int base_size, const int* pair_connection, int pair_size, int* int_seq)
{
   double energy = 0.0;
   double min;
   int i;

   for (i=0; i<base_size; i++)
   {
      min = MAX_DOUBLE;

      /*There are 4 cases:  (a) BP b BP b
                            (b) b BP b BP b
                            (c) b BP b BP
                            (d) BP b BP*/

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         // corresponds to cases (b) and (c) (from the second free base on)
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i-1]][1]]+4*int_seq[BP_Order[pair_connection[i-1]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i-1]][1]]+4*int_seq[BP_Order[pair_connection[i-1]][0]]+int_seq[base_connection[i]]];
      }

      if (i < pair_size)
      {
         // usual case: base is located left of (= before) the stem
         // corresponds to cases (b) and (c) (free bases in front of a following stem)
         if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
            if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]];

         // usual case: base is located right of (= after) the stem
         // corresponds to cases (a) and (d) (free bases after the stem)
         if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
            if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]];
      }

      if ( i+1 < pair_size)
      {
         // usual case: base is located left of (=before) the stem
         // correponds to cases (a) and (d) (free base in front of the stem)
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i+1]][1]]+4*int_seq[BP_Order[pair_connection[i+1]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i+1]][1]]+4*int_seq[BP_Order[pair_connection[i+1]][0]]+int_seq[base_connection[i]]];
      }
      energy = Sum_MaxDouble3(energy, min, BasePenalty(base_connection[i], int_seq[base_connection[i]]));
   }
   return energy;
}



/****************************************************************************************

calculates the MINIMAL (=best) free energy of the external loop, values are stored
in the last row of D and Trace

****************************************************************************************/

void externBestEnergy()
{
   // in the external loop we have to take into account the ending stems, the free base between them, and the dangling ends
   /*e.g.:    bp_i - bp_j                    bp_i - bp_j
         closing_i - closing_j freebase closing_i - closing_j free_base*/

   //finding all predecessors of the external loop:
   /*checking the structure vector (brackets), start at the first opening bracket, jump to its closing counterpart (or to the following pos.),
     number of stems in the external loop++, then looking for the next opening bracket, jump to its closing counterpart +1, ....*/
   int bp_i, bp_j;
   unsigned int i;

   int* vorg_last; // BP_Order-pos. of the BPs that are predecessors of the external loop
   vorg_last = (int*) malloc(sizeof(int)*BP_Order[numBP][3]);
   for (int a=0; a<BP_Order[numBP][3]; a++)
      vorg_last[a] = -1;

   int stems = 0;
   int pos = 0;
   while (pos < struct_len)
   {
      if (brackets[pos] == '(')
      {
         vorg_last[stems] = BP_Pos_Nr[pos];
         stems++;
         pos = BP_Order[BP_Pos_Nr[pos]][1] + 1;
      }
      else
         pos++;
   }

   //initially, the last row of D is filled with values of the row before, if there are no dang-ends, no further changes have to be done to the last row
   for (int bp_assign = 0; bp_assign < 6; bp_assign++)
      D[numBP][bp_assign] = D[numBP-1][bp_assign];

   double min_dang;
   double front_dang = 0.0;
   double base_energy;

   vector< vector<int> > EndBaseConnections;     // all connections, concerning free bases
   vector< vector<int> > EndBasePairConnections; // all connections, concerning stems (BP)
   EndConnections(Maximum(BP_Order[numBP][3],1), vorg_last, EndBaseConnections, EndBasePairConnections);

   for (int bp_assign = 0; bp_assign<6; bp_assign++)
   {
      for (int vg = 0; vg < BP_Order[numBP][3]; vg++)
         Trace[numBP][bp_assign][vg][0] = vorg_last[vg];
      //during the traceback go to the same assignment as the position in the artificial last row of D and Trace 
      Trace[numBP][bp_assign][0][1] = bp_assign; //--> this can be overwritten later, but if the second if-loop is not reached, this initialization is needed

      // if there are free bases in front of the first opening bracket
      if (BP_Order[numBP-1][0] > 0)
      {
         min_dang = MAX_DOUBLE;
         BP2_2(bp_assign, bp_i, bp_j);
         for (int free_base = 0; free_base < 4; free_base++)
         {
            base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_j+4*bp_i+free_base], BasePenalty(BP_Order[numBP-1][0]-1,free_base));
            if (min_dang > base_energy)
               min_dang = base_energy;
         }
         front_dang = min_dang;
      }

      // if there is a structural fragment after the closing bracket of the last BP:
      if (BP_Order[numBP-1][1] < struct_len-1)
      {
         double min = MAX_DOUBLE;
         double energy=0.0;
         double energy_help = 0.0, energy_unbound, energy_stems;
         int pair_size, base_size;

         int** bp_at_pair_connection;  // base pair assignment for a connection
         int** MIN_bp_at_pair_connection;  // base pair assignment for a connection with min. free energy

         for (i=0; i<EndBaseConnections.size(); i++) //EndBaseConnections.size() == EndBasePairConnections.size()
         {
            pair_size = EndBasePairConnections[i].size();
            base_size = EndBaseConnections[i].size();

            bp_at_pair_connection = (int**) malloc(sizeof(int*)*pair_size);
            MIN_bp_at_pair_connection = (int**) malloc(sizeof(int*)*pair_size);
            for (int a=0; a<pair_size; a++)
            {
               bp_at_pair_connection[a] = (int*) malloc(sizeof(int)*2);
               MIN_bp_at_pair_connection[a] = (int*) malloc(sizeof(int)*2);
               for (int b=0; b<2; b++)
               {
                   bp_at_pair_connection[a][b] = 0;
                   MIN_bp_at_pair_connection[a][b] = 0;
               }
            }

            /* for bp_at_pair_connection, each possible assignment has to be tested. Thereto a vector of size 'pair_size' is generated

            the assignments are generated by generating an integer value that count up to the number of possible assignments (6^pair_size)
            and then evaluated the representation in the 6-Number-System, e.g.: 345 = 0*6^4+1*6^3+3*6^2+3*6^1+3*6^0 ^= 01333*/

            // before, it is tested whether the last BP (conc. BP_Order) is included in this connection. Since this is already fixed,
            // no further variation is needed (if the last BP is included in the connection, then at the first pos.)

            int start_pos; // if the last BP is the first BP in the connection, then the start_pos for the assignment of bp_at_pair_connection
                           // is not 0 but 1
            int start_power;

            if (EndBasePairConnections[i][0] == numBP-1)
            {
               start_pos = 1;
               BP2_2(bp_assign, bp_i, bp_j);
               bp_at_pair_connection[0][0] = bp_i;
               bp_at_pair_connection[0][1] = bp_j;
               start_power = pair_size-2;
            }
            else
            {
               start_pos = 0;
               start_power = pair_size-1;
            }

            min = MAX_DOUBLE;

            for (unsigned long int int_anz = 0; int_anz < (unsigned long int)pow(6,(double)(start_power+1)); int_anz++)
            {
               int power = start_power;
               int value = int_anz;
               int pos = start_pos;
               int bp_help;

               // create an assignment
               while (power >= 0)
               {
                  bp_help = value/(int)pow(6,(double)power);

                  BP2_2(bp_help, bp_i, bp_j);
                  bp_at_pair_connection[pos][0]=bp_i;
                  bp_at_pair_connection[pos][1]=bp_j;

                  value = (int)fmod((double)value,pow(6,(double)power));

                  power--;
                  pos++;
               }

               // energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here, energy of 
               // the stem will be added later
               energy_unbound = BestEndConnectionEnergy(EndBaseConnections[i], EndBasePairConnections[i], (const int**) bp_at_pair_connection);

               // add the energy of the stems (depending on the the closing BP), BUT: not closing BP!
               energy_stems = 0.0;
               for (int j=0; j<pair_size; j++)
                  energy_stems = Sum_MaxDouble(energy_stems,D[EndBasePairConnections[i][j]][BP2int(bp_at_pair_connection[j][0],bp_at_pair_connection[j][1])]);

               energy_help = Sum_MaxDouble(energy_stems,energy_unbound);

               if (min > energy_help)
               {
                  min = energy_help;
                  for (int a=0; a<pair_size; a++)
                  {
                      MIN_bp_at_pair_connection[a][0] = bp_at_pair_connection[a][0];
                      MIN_bp_at_pair_connection[a][1] = bp_at_pair_connection[a][1];
                  }
               }
            }// for int_anz

            // if no valid assignment was found, set a random one (this is the case, if in the external loop a BP exists, that has no valid assignment
            // concerning the constraints --> usually this error is filtered before)
            if (min == MAX_DOUBLE)
               for (int a=0; a<pair_size; a++)
               {
                  int rand_bp = RandomBasePair();
                  BP2_2(rand_bp, MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);
               }

            energy = Sum_MaxDouble(energy,min);

            // store best assignment of the predecessors in the traceback
            for (int a=0; a<pair_size; a++)
               for (int vg = 0; vg<BP_Order[numBP][3]; vg++)
                  if (EndBasePairConnections[i][a] == Trace[numBP][bp_assign][vg][0])
                     Trace[numBP][bp_assign][vg][1] = BP2int(MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);

            for (int c=0; c<pair_size; c++)
            {
               free(bp_at_pair_connection[c]);
               free(MIN_bp_at_pair_connection[c]);
            }
            free(bp_at_pair_connection);
            free(MIN_bp_at_pair_connection);
         } // for i

         D[numBP][bp_assign] = Sum_MaxDouble(energy,front_dang);
      //}// else
      }//if

   } //for bp_assign

   for (i=0; i<EndBaseConnections.size(); i++)
      EndBaseConnections[i].clear();
   EndBaseConnections.clear();

   for (i=0; i<EndBasePairConnections.size(); i++)
      EndBasePairConnections[i].clear();
   EndBasePairConnections.clear();
}


/****************************************************************************************

 finds energy of the external loop (for a given sequence!!)

****************************************************************************************/

double externEnergy(int* int_seq)
{
   // in the external loop we have to take into account the ending stems, the free base between them, and the dangling ends
   /*e.g.:    bp_i - bp_j                    bp_i - bp_j
         closing_i - closing_j freebase closing_i - closing_j free_base*/

   int bp_i, bp_j, a, free_base;
   int i, j;

   int stems = 0;
   int pos = 0;

   int** EndBaseConnections;     // all connections, concerning free bases
   int** EndBasePairConnections; // all connections, concerning stems
   int* EndBaseSizes;            // sizes of the base connections
   int* EndBasePairSizes;        // sizes of the stem connections
   int* vorg_last;               // BP_Order-pos. of the BPs that are predecessors of the external loop

   int num_conns;
   double energy=0.0;
   double energy_unbound;
   int pair_size, base_size;

   int con_size = BP_Order[numBP][3]+2;  // in order that base and basepairconn. have the same size, no conn. can be bigger than #stems+2 and all together
                                         // there can not be more (since only adjacent free base are taken into account)

   //finding all predecessors of the external loop:
   /*checking the structure vector (brackets), start at the first opening bracket, jump to its closing counterpart (or to the following pos.),
     number of stems in the external loop++, then looking for the next opening bracket, jump to its closing counterpart +1, ....*/
   vorg_last = (int*) malloc(sizeof(int)*BP_Order[numBP][3]);
   for (a=0; a<BP_Order[numBP][3]; a++)
      vorg_last[a] = -1;

   while (pos < struct_len)
   {
      if (brackets[pos] == '(')
      {
         vorg_last[stems++] = BP_Pos_Nr[pos];
         pos = BP_Order[BP_Pos_Nr[pos]][1] + 1;
      }
      else
         pos++;
   }

   // if there are free bases in front of the first opening bracket
   if (BP_Order[numBP-1][0] > 0)
   {
      bp_i = int_seq[BP_Order[numBP-1][0]];
      bp_j = int_seq[BP_Order[numBP-1][1]];
      free_base = int_seq[BP_Order[numBP-1][0]-1];
      energy = Sum_MaxDouble3(energy, single_base_stacking_energy[64+16*bp_j+4*bp_i+free_base], BasePenalty(BP_Order[numBP-1][0]-1,free_base));
      //if there are more free bases in front of the frist opening bracket, the penalties are considered seperately in get_BasePart_Energy in 
      //multi_energy.cpp
   }

   // since these arrays here can not be created dynamically, an array of maximal size is created
   EndBaseConnections = (int**) malloc(sizeof(int*)*con_size);
   EndBaseSizes = (int*) malloc(sizeof(int)*con_size);
   for (i=0; i<con_size; i++)
   {
      EndBaseConnections[i] = (int*) malloc(sizeof(int) * con_size);
      EndBaseSizes[i] = -1;
      for (j=0;j<con_size; j++)
         EndBaseConnections[i][j] = -1;
   }

   EndBasePairConnections = (int**) malloc(sizeof(int*)*con_size);
   EndBasePairSizes = (int*) malloc(sizeof(int) * con_size);
   for (i=0; i<con_size; i++)
   {
      EndBasePairConnections[i] = (int*) malloc(sizeof(int) * con_size);
      EndBasePairSizes[i] = -1;
      for (j=0;j<con_size; j++)
         EndBasePairConnections[i][j] = -1;
   }

   EndConnections(vorg_last, EndBaseConnections, EndBaseSizes, EndBasePairConnections, EndBasePairSizes);

   //count number of connections
   num_conns = 0;
   for (i=0; i<Maximum(BP_Order[numBP][2],BP_Order[numBP][3]+2); i++)
   {
      if ((EndBaseConnections[i][0] != -1) || (EndBasePairConnections[i][0] != -1))
         num_conns++;
      else
         break;
   }

         for (i=0; i<num_conns; i++) //EndBaseConnections.size() == EndBasePairConnections.size()
         {
            pair_size = EndBasePairSizes[i];
            base_size = EndBaseSizes[i];

            // energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here, energy of 
            // the stem will be added later
            energy_unbound = EndConnectionEnergy(EndBaseConnections[i], base_size, EndBasePairConnections[i], pair_size, int_seq);
            energy = Sum_MaxDouble(energy, energy_unbound);
         }
   for (i=0; i<con_size; i++)
      free(EndBaseConnections[i]);
   free(EndBaseConnections);

   for (i=0; i<con_size; i++)
      free(EndBasePairConnections[i]);
   free(EndBasePairConnections);

   return energy;
}



/******************************************************

finds the assignment with minimal free energy for the 
free bases for one connection if the assignments of the
last BPs of the stems are known
(only needed for Traceback)

base_connection = one connection (conc. free bases)
pair_connection = one connection (conc. stems)
bp_at_pair_connection = BP assignment for pair_connection

******************************************************/

int* EndConnectionBestFreeBases(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection)
{
   double min = MAX_DOUBLE;
   double base_energy;
   int base_size = (int)base_connection.size();
   int pair_size = (int)pair_connection.size();
   int min_base=0;
   int* min_bases;
   min_bases = (int*) malloc(sizeof(int)*base_size);
   for (int i=0; i<base_size; i++)
      min_bases[i] = -1;

   for (int i=0; i<base_size; i++)
   {
      //for each base: testing whether it is adjacent to a stem or two, minimize the energy
      min = MAX_DOUBLE;

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }
      }

      if (i < pair_size)
      {
         // usual case: base is located left of (= before) the stem
         if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }

         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }
      }

      if ( i+1 < (int)pair_connection.size())
      {
         // usual case: base is located left of (=before) the stem
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][1]+4*bp_at_pair_connection[i+1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }
      }
      min_bases[i] = min_base;
   }
   return min_bases;
}






