#include "multi_energy.h"


/*******************************************************
 identifies connections concerning the closing base pairs
 of the stems and concerning the free bases adjacent to them

 stem_ends = contains all stem ends and their positions in BP_Order, we are interested in the last stem_num
 order_pos_of_closing_bp = pos. of the closing BP in BP_Order
 BaseConnections = contains row-wise the connections concerning bases
 BasePairConnections = contains row-wise the connections concerning base pairs

 Return: vector including BPs in the ML (= predecessors of the closing BP of the ML)
*******************************************************/

int* MultiLoopConnections(int stem_num, vector<int> & stem_ends, int order_pos_of_closing_bp, vector< vector<int> > & BaseConnections, vector< vector<int> > & BasePairConnections)
{
   int* ml_stem_ends = (int*) malloc(sizeof(int)*stem_num);      //stem ends in the ML
   int* ml_between_stems = (int*) malloc(sizeof(int)*stem_num+1);  //contains the number of bases between the stems (in each coord. the bases "before" the
                                                                    //stem, only at the last coord. of bases_between_stems the free bases after
                                                                    //the last stem and before the closing BP are stored

   vector<int> one_connection;         //help vector for one connection

   // get just the interesting stem ends from stem_ends and erase them there,
   // thereby their order is reversed. BPs (i,j) having a small i are first now
   for (int i=0; i<stem_num; i++)
   {
      ml_stem_ends[i] = stem_ends.back();
      stem_ends.pop_back();
   }

   // number of free bases between the stem-ending base pairs
   ml_between_stems[0] = BP_Order[ml_stem_ends[0]][0] - BP_Order[order_pos_of_closing_bp][0] - 1;
   for (int i=1; i<stem_num; i++)
      ml_between_stems[i] = BP_Order[ml_stem_ends[i]][0] - BP_Order[ml_stem_ends[i-1]][1] - 1;
   ml_between_stems[stem_num] = BP_Order[order_pos_of_closing_bp][1] - BP_Order[ml_stem_ends[stem_num-1]][1] - 1;

   // finding all connections: (concerning free bases):
   //*****************************************************

   // if there's at least one free base upstream the first stem
   if (ml_between_stems[0] != 0)
      one_connection.push_back(BP_Order[order_pos_of_closing_bp][0]+1);

   for (int i=0; i<stem_num; i++)
   {
      // if there is only one free base between the stems and if there is at least one free base after the next stem
      if (ml_between_stems[i] == 1)
      {
         if (ml_between_stems[i+1] != 0)
            one_connection.push_back(BP_Order[ml_stem_ends[i]][1]+1);
      }
      else
      {
         BaseConnections.push_back(one_connection);
         one_connection.clear();
         if (ml_between_stems[i] != 0)
            one_connection.push_back(BP_Order[ml_stem_ends[i]][0]-1);
         if (ml_between_stems[i+1] != 0)
            one_connection.push_back(BP_Order[ml_stem_ends[i]][1]+1);
      }
   }

   if (ml_between_stems[stem_num] == 1)
   {
      BaseConnections.push_back(one_connection);
      one_connection.clear();
   }
   else
   {
      BaseConnections.push_back(one_connection);
      one_connection.clear();
      if (ml_between_stems[stem_num] != 0)
      {
         one_connection.push_back(BP_Order[order_pos_of_closing_bp][1]-1);
      }
      // one_connection might be empty, but this is done to ensure that BaseConnections and BasePairConnections
      // are always of the same size
      BaseConnections.push_back(one_connection);
      one_connection.clear();
   }

   // finding all connections: (concerning stems):
   //***********************************************

   one_connection.push_back(order_pos_of_closing_bp);

   for (int i=0; i<stem_num; i++)
   {
      if (ml_between_stems[i] == 1)
         one_connection.push_back(ml_stem_ends[i]);
      else
      {
         BasePairConnections.push_back(one_connection);
         one_connection.clear();
         one_connection.push_back(ml_stem_ends[i]);
      }
   }

   if (ml_between_stems[stem_num] == 1)
   {
      one_connection.push_back(order_pos_of_closing_bp);
      BasePairConnections.push_back(one_connection);
      one_connection.clear();
   }
   else
   {
      BasePairConnections.push_back(one_connection);
      one_connection.clear();
      one_connection.push_back(order_pos_of_closing_bp);
      BasePairConnections.push_back(one_connection);
      one_connection.clear();
   }

   return ml_stem_ends;
}




/*******************************************************
 identifies connections concerning the closing base pairs
 of the stems and concerning the free bases adjacent to them
 (similar to the other MultiLoopConnections-function, but without
 using C++ vectors)

 stem_ends = contains all stem ends of leaving stems and their positions in BP_Order
 order_pos_of_closing_bp = pos. of the closing BP in BP_Order
 BaseConnections = contains row-wise the connections concerning bases
 BasePairConnections = contains row-wise the connections concerning base pairs
 EndBaseSizes = sizes of the EndBaseConnections
 EndBasePairSizes = sizes of the EndBasePairConnections
*******************************************************/

void MultiLoopConnections(int* stem_ends, int order_pos_of_closing_bp, int** BaseConnections, int *BaseSizes, int** BasePairConnections, int* BasePairSizes)
{
   int i;
   int stem_num = BP_Order[order_pos_of_closing_bp][3];
   int num_bp_con = 0;   //current number of BasePairConnections
   int num_base_con = 0; //current number of BaseConnections
   int bp, base;         //current position in the current connection

   int* ml_stem_ends = (int*) malloc(sizeof(int)*stem_num);  //stem ends of the ML
   int* ml_between_stems = (int*) malloc(sizeof(int)*stem_num+1);  //contains the number of bases between the stems (in each coord. the bases "before" the
                                                                   //stem, only at the last coord. of bases_between_stems the free bases after
                                                                   //the last stem and before the closing BP are stored

   // get just the interesting stem ends from stem_ends and erase them there,
   // thereby their order is reversed. BPs (i,j) having a small i are first now
   for (i=0; i<stem_num; i++)
      ml_stem_ends[i] = stem_ends[i];

   // number of free bases between the stem-ending base pairs
   ml_between_stems[0] = BP_Order[ml_stem_ends[0]][0] - BP_Order[order_pos_of_closing_bp][0] - 1;
   for (i=1; i<stem_num; i++)
      ml_between_stems[i] = BP_Order[ml_stem_ends[i]][0] - BP_Order[ml_stem_ends[i-1]][1] - 1;
   ml_between_stems[stem_num] = BP_Order[order_pos_of_closing_bp][1] - BP_Order[ml_stem_ends[stem_num-1]][1] - 1;


   // finding all connections: (concerning free bases):
   //*****************************************************

   bp = 0;
   base = 0;

   // if there's at least one free base upstream the first stem
   if (ml_between_stems[0] != 0)
      BaseConnections[num_base_con][base++] = BP_Order[order_pos_of_closing_bp][0]+1;


   for (i=0; i<stem_num; i++)
   {
      // if there is only one free base between the stems and if there is at least one free base after the next stem
      if (ml_between_stems[i] == 1)
      {
         if (ml_between_stems[i+1] != 0)
            BaseConnections[num_base_con][base++] = BP_Order[ml_stem_ends[i]][1]+1;
      }
      else
      {
         BaseSizes[num_base_con] = base;
         num_base_con++;
         base = 0;
         if (ml_between_stems[i] != 0)
            BaseConnections[num_base_con][base++] = BP_Order[ml_stem_ends[i]][0]-1;
         if (ml_between_stems[i+1] != 0)
            BaseConnections[num_base_con][base++] = BP_Order[ml_stem_ends[i]][1]+1;
      }
   }

   if (ml_between_stems[stem_num] == 1)
   {
      BaseSizes[num_base_con] = base;
      num_base_con++;
      base = 0;
   }
   else
   {
      BaseSizes[num_base_con] = base;
      num_base_con++;
      base = 0;
      if (ml_between_stems[stem_num] != 0)
         BaseConnections[num_base_con][base++] = BP_Order[order_pos_of_closing_bp][1]-1;
      // one_connection might be empty, but this is done to ensure that BaseConnections and BasePairConnections
      // are always of the same size
      BaseSizes[num_base_con] = base;
      num_base_con++;
      base = 0;
   }

   // finding all connections: (concerning stems):
   //************************************************

   BasePairConnections[num_bp_con][bp++] = order_pos_of_closing_bp;

   for (i=0; i<stem_num; i++)
   {
      if (ml_between_stems[i] == 1)
         BasePairConnections[num_bp_con][bp++] = ml_stem_ends[i];
      else
      {
         BasePairSizes[num_bp_con] = bp;
         num_bp_con++;
         bp = 0;
         BasePairConnections[num_bp_con][bp++] = ml_stem_ends[i];
      }
   }

   if (ml_between_stems[stem_num] == 1)
   {
      BasePairConnections[num_bp_con][bp++] = order_pos_of_closing_bp;
      BasePairSizes[num_bp_con] = bp;
      num_bp_con++;
      bp = 0;
   }
   else
   {
      BasePairSizes[num_bp_con] = bp;
      num_bp_con++;
      bp = 0;
      BasePairConnections[num_bp_con][bp++] = order_pos_of_closing_bp;
      BasePairSizes[num_bp_con] = bp;
      num_bp_con++;
      bp = 0;
   }
}




/******************************************************
gives the MINIMAL (=best) free energy for a connection of the
ML if the base pair assignments of the stem-ending
base pairs are known

base_connection = a connection concerning free bases
pair_connection = a connection concerning stem-ending base pairs
bp_at_pair_connection = a BP assignment of pair_connection
******************************************************/

double BestConnectionEnergy(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection)
{
   double energy = 0.0;
   double min = MAX_DOUBLE;
   double base_energy;

   for (int i=0; i<(int)base_connection.size(); i++)
   {
      //for each base, test, whether it is adjacent to one or two stems, then: minimize energy

      //bases in base_connection[i] may be neighbored to base pairs in
      //pair_connection[i-1] [i] or/and [i+1] => testing (and also whether the respecting field in the vector exists)
      min = MAX_DOUBLE;

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }
      }

      if (i < (int)pair_connection.size())
      {
         // if the left stem is the closingBP
         if (BP_Order[pair_connection[i]][0]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located left of (= before) the stem
         if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // if the right stem is the closingBP
         if (BP_Order[pair_connection[i]][1]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

      }

      if ( i+1 < (int)pair_connection.size())
      {
         // usual case: base is located left of (= before) the stem
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][1]+4*bp_at_pair_connection[i+1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i+1]][1]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][0]+4*bp_at_pair_connection[i+1][1]+j],BasePenalty(base_connection[i],j));
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

double ConnectionEnergy(const int* base_connection, int base_size, const int* pair_connection, int bp_size, int* int_seq)
{
   double energy = 0.0;
   int i;
   double min;

   for (i=0; i<base_size; i++)
   {
      //for each base, test, whether it is adjacent to one or two stems, then: minimize energy

      //bases in base_connection[i] may be neighbored to base pairs in
      //pair_connection[i-1] [i] or/and [i+1] => testing (and also whether the respecting field in the vector exists)
      min = MAX_DOUBLE;

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i-1]][1]]+4*int_seq[BP_Order[pair_connection[i-1]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i-1]][1]]+4*int_seq[BP_Order[pair_connection[i-1]][0]]+int_seq[base_connection[i]]];
      }

      if (i < bp_size)
      {
         // if the left stem is the closingBP
         if (BP_Order[pair_connection[i]][0]+1 == base_connection[i])
            if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][0]]+4*int_seq[BP_Order[pair_connection[i]][1]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][0]]+4*int_seq[BP_Order[pair_connection[i]][1]]+int_seq[base_connection[i]]];

         // usual case: base is located left of (=before) the stem
         if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
            if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]];

         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
            if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]];

         // if the right stem is the closingBP
         if (BP_Order[pair_connection[i]][1]-1 == base_connection[i])
            if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][0]]+4*int_seq[BP_Order[pair_connection[i]][1]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][0]]+4*int_seq[BP_Order[pair_connection[i]][1]]+int_seq[base_connection[i]]];
      }

      if ( i+1 < bp_size)
      {
         // usual case: base is located left of (=before) the stem
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i+1]][1]]+4*int_seq[BP_Order[pair_connection[i+1]][0]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i+1]][1]]+4*int_seq[BP_Order[pair_connection[i+1]][0]]+int_seq[base_connection[i]]];

         // if the right stem is the closingBP
         if (BP_Order[pair_connection[i+1]][1]-1 == base_connection[i])
            if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i+1]][0]]+4*int_seq[BP_Order[pair_connection[i+1]][1]]+int_seq[base_connection[i]]] < min)
               min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i+1]][0]]+4*int_seq[BP_Order[pair_connection[i+1]][1]]+int_seq[base_connection[i]]];
      }
      //base penalty
      energy = Sum_MaxDouble3(energy, min, BasePenalty(base_connection[i],int_seq[base_connection[i]]));
   }
   return energy;
}





/****************************************************************************************

calculates the MINIMAL (=best) free energy of the ML, values are stored in D and Trace

****************************************************************************************/

void MLBestEnergy(int bp_pos, vector<int> & stem_ends)
{
   double min = MAX_DOUBLE;
   double MLenergy, energy;
   double offset = 3.4;
   double free_base_penalty = 0.0;
   double helix_penalty = 0.4;
   double energy_help = 0.0, energy_unbound, energy_stems;
   int bp_i, bp_j;
   int pair_size, base_size;
   unsigned int i;
   int* ML_vorgaenger; // BP_Order pos. of the stems in the ML
   vector< vector<int> > BaseConnections;     // all connections conc. free bases
   vector< vector<int> > BasePairConnections; // all connections conc. stems

   int** bp_at_pair_connection;      // base pair assignment for a connection
   int** MIN_bp_at_pair_connection;  // base pair assignment for a connection with min. free energy

   stem_ends.push_back(bp_pos-1);    //previous BP is a stem ending BP in the ML
   //BP_Order[bp_pos][3]+1, since also closingBP-stems counts
   MLenergy = offset + free_base_penalty * BP_Order[bp_pos][2] + helix_penalty * (BP_Order[bp_pos][3]+1);
   ML_vorgaenger = MultiLoopConnections(BP_Order[bp_pos][3], stem_ends, bp_pos, BaseConnections, BasePairConnections);

   // since the closing BP of the ML can be included in two connections, its assignment is fixed before testing the assignments of the other BPs
   for (int bp_assign=0; bp_assign<6; bp_assign++)
   {
      //store predecessors in Trace, assignment still unknown
      for (int vg=0; vg<BP_Order[bp_pos][3]; vg++)
         Trace[bp_pos][bp_assign][vg][0] = ML_vorgaenger[vg];

      energy = MLenergy;  // energy

      for (i=0; i<BaseConnections.size(); i++) //BaseConnections.size() == BasePairConnections.size()
      {
         pair_size = BasePairConnections[i].size();
         base_size = BaseConnections[i].size();

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

         int start_pos; // if the closing BP is the first one in the connection, start_pos for assigning bp_at_pair_connection
                        // is not 0, but 1
         int start_power;  // if the closing BP is the last in the connection, start_power for assigning bp_at_pair_connection
                           // is not size-1, but size-2, since then the last position is already fixed
         if (BasePairConnections[i][0] == bp_pos)
         {
            start_pos = 1;
            BP2_2(bp_assign, bp_i, bp_j);
            bp_at_pair_connection[0][0] = bp_i;
            bp_at_pair_connection[0][1] = bp_j;

            if (BasePairConnections[i][pair_size-1] == bp_pos)
            {
               start_power = pair_size-3;
               bp_at_pair_connection[pair_size-1][0] = bp_i;
               bp_at_pair_connection[pair_size-1][1] = bp_j;
            }
            else
               start_power = pair_size-2;
         }
         else
         {
            start_pos = 0;
            if (BasePairConnections[i][pair_size-1] == bp_pos)
            {
               start_power = pair_size-2;
               BP2_2(bp_assign, bp_i, bp_j);
               bp_at_pair_connection[pair_size-1][0] = bp_i;
               bp_at_pair_connection[pair_size-1][1] = bp_j;
            }
            else
               start_power = pair_size-1;
         }

         min = MAX_DOUBLE;

         //this is the case, if the connection only consists of the closing BP, then the following for-loop is not done.
         //Nevertheless, an adjacent free base has to be analyzed.
         if ((pair_size == 1) && (start_power == pair_size-3))
         {
            if (base_size == 0)
               min = 0;
            else
            {
               // energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here)
               energy_unbound = BestConnectionEnergy(BaseConnections[i], BasePairConnections[i], (const int**) bp_at_pair_connection);
               min = energy_unbound;
            }
         }

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
            energy_unbound = BestConnectionEnergy(BaseConnections[i], BasePairConnections[i], (const int**) bp_at_pair_connection);

            // add the energy of the stems (depending on the the closing BP), BUT: not closing BP!
            energy_stems = 0.0;
            for (int j=0; j<pair_size; j++)
               if (BasePairConnections[i][j] != bp_pos)
                  energy_stems = Sum_MaxDouble(energy_stems,D[BasePairConnections[i][j]][BP2int(bp_at_pair_connection[j][0],bp_at_pair_connection[j][1])]);

            energy_help = Sum_MaxDouble(energy_stems, energy_unbound);

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

         // if no valid assignment was found, set a random one (this is the case, if in the loop a BP exists, that has no valid assignment
         // concerning the constraints --> usually this error is filtered before)
         if (min == MAX_DOUBLE)
            for (int a=0; a<pair_size; a++)
            {
               int rand_bp = RandomBasePair();
               BP2_2(rand_bp, MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);
            }

         // add minimal energy that arises because of the free bases in the ML
         energy = Sum_MaxDouble(energy,min);

         // store best assignment of the predecessors in the traceback
         for (int a=0; a<pair_size; a++)
            if (BasePairConnections[i][a] != bp_pos)
               for (int vg = 0; vg<BP_Order[bp_pos][3]; vg++)
                  if (BasePairConnections[i][a] == Trace[bp_pos][bp_assign][vg][0])
                     Trace[bp_pos][bp_assign][vg][1] = BP2int(MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);

         for (int c=0; c<pair_size; c++)
         {
            free(bp_at_pair_connection[c]);
            free(MIN_bp_at_pair_connection[c]);
         }
         free(bp_at_pair_connection);
         free(MIN_bp_at_pair_connection);
      } // for i
      
      BP2_2(bp_assign, bp_i, bp_j);
      D[bp_pos][bp_assign] = Sum_MaxDouble3(energy, Zero_or_StemEndAU(bp_pos, bp_assign),PairPenalty(bp_pos, bp_i, bp_j));

   } //for bp_assign

   for (i=0; i<BaseConnections.size(); i++)
      BaseConnections[i].clear();
   BaseConnections.clear();

   for (i=0; i<BasePairConnections.size(); i++)
      BasePairConnections[i].clear();
   BasePairConnections.clear();
}




/****************************************************************************************

finds energy of the ML (for a given sequence!!) (only the fraction that is relevant
for the difference of two assignments of it)

****************************************************************************************/

double MLEnergy(int bp_pos, int* int_seq)
{
   double energy = 0.0;
   //double offset = 3.4;
   //double free_base_penalty = 0.0;
   //double helix_penalty = 0.4;
   double energy_unbound;
   int pair_size, base_size;
   unsigned int i, j;
   int** BaseConnections;     // all connections, concerning free bases
   int** BasePairConnections; // all connections, concerning stems
   int* BaseSizes;            // sizes of the base connections
   int* BasePairSizes;        // sizes of the stem connections

   int bp_assign; // assignment of the closing BP
   int num_conns;
   int* stem_ends;
   unsigned int con_size = BP_Order[bp_pos][3]+2; // in order that base and basepairconn. have the same size, no conn. can be bigger than #stems+2 and all together
                                                  // there can not be more (since only adjacent free base are taken into account)
   stem_ends = FindStemEnds(bp_pos);

   BaseConnections = (int**) malloc(sizeof(int*)*con_size);
   BaseSizes = (int*) malloc(sizeof(int)*con_size);
   for (i=0; i<con_size; i++)
   {
      BaseConnections[i] = (int*) malloc(sizeof(int) * con_size);
      BaseSizes[i] = -1;
      for (j=0;j<con_size; j++)
         BaseConnections[i][j] = -1;
   }

   BasePairConnections = (int**) malloc(sizeof(int*)*con_size);
   BasePairSizes = (int*) malloc(sizeof(int) * con_size);
   for (i=0; i<con_size; i++)
   {
      BasePairConnections[i] = (int*) malloc(sizeof(int) * con_size);
      BasePairSizes[i] = -1;
      for (j=0;j<con_size; j++)
         BasePairConnections[i][j] = -1;
   }

   //BP_Order[bp_pos][3]+1, since also the closing BP counts
   //energy = offset + free_base_penalty * BP_Order[bp_pos][2] + helix_penalty * (BP_Order[bp_pos][3]+1);
   MultiLoopConnections(stem_ends, bp_pos, BaseConnections, BaseSizes, BasePairConnections, BasePairSizes);
   bp_assign = BP2int(int_seq[BP_Order[bp_pos][0]], int_seq[BP_Order[bp_pos][1]]);

   //count the existing connections:
   num_conns = 0;
   for (i=0; i<con_size; i++)
   {
      if ((BaseConnections[i][0] != -1) || (BasePairConnections[i][0] != -1))
         num_conns++;
      else
         break;
   }

   for (i=0; i<(unsigned int)num_conns; i++) //BaseConnections.size() == BasePairConnections.size()
   {
      pair_size = BasePairSizes[i];
      base_size = BaseSizes[i];

      // energy of one connection (just energy of the free bases)
      energy_unbound = ConnectionEnergy(BaseConnections[i], base_size, BasePairConnections[i], pair_size, int_seq);
      energy = Sum_MaxDouble(energy, energy_unbound);
   } // for i

   //base pair penalty only for the closingBP
   energy = Sum_MaxDouble3(energy, PairPenalty(bp_pos,int_seq[BP_Order[bp_pos][0]],int_seq[BP_Order[bp_pos][1]]), Zero_or_StemEndAU(bp_pos, bp_assign));

   for (i=0; i<con_size; i++)
      free(BaseConnections[i]);
   free(BaseConnections);

   for (i=0; i<con_size; i++)
      free(BasePairConnections[i]);
   free(BasePairConnections);

   return energy;
}



/*****************************************************************

Finds all BPs that are stems ends in a ML

*****************************************************************/

int* FindStemEnds(int closingBP)
{
   int i, bp_j;
   int* stem_ends;
   stem_ends = (int*) malloc(sizeof(int)*BP_Order[closingBP][3]);

   //previous base pair in BP_Order is a closingBP
   stem_ends[0] = closingBP - 1;

   for (i=1; i<BP_Order[closingBP][3]; i++)
   {
      bp_j = BP_Order[stem_ends[i-1]][1]+1;
      while (brackets[bp_j] != '(')
         bp_j++;
      stem_ends[i] = BP_Pos_Nr[bp_j];
   }
   return stem_ends;
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

int* ConnectionBestFreeBases(const vector<int> & base_connection, const vector<int> & pair_connection, const int** bp_at_pair_connection)
{
   double min = MAX_DOUBLE;
   double base_energy;
   int base_size = (int)base_connection.size();
   int min_base=0;
   int* min_bases;
   min_bases = (int*) malloc(sizeof(int)*base_size);
   for (int i=0; i<base_size; i++)
      min_bases[i] = -1;

   for (int i=0; i<base_size; i++)
   {
      //for each base, test, whether it is adjacent to one or two stems, then: minimize energy

      //bases in base_connection[i] may be neighbored to base pairs in
      //pair_connection[i-1] [i] or/and [i+1] => testing (and also whether the respecting field in the vector exists)
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

      if (i < (int)pair_connection.size())
      {
         // if the left stem is the closing BP
         if (BP_Order[pair_connection[i]][0]+1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }

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

         // if the right stem is the closing BP
         if (BP_Order[pair_connection[i]][1]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }
      }

      if ( i+1 < (int)pair_connection.size())
      {
         // usual case: base is located left of (= before) the stem
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

         // if the right stem is the closing BP
         if (BP_Order[pair_connection[i+1]][1]-1 == base_connection[i])
            for (int j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][0]+4*bp_at_pair_connection[i+1][1]+j],BasePenalty(base_connection[i],j));
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



