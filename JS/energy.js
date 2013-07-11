/***********************************************************
README:

MultiLoopConnections is overridden in the c++ code, something we can't do in javascript, so
I made 2 seperate functions : MultiLoopConnections (The first one in multi_energy.cpp) and 
MultiLoopConnections2 (the second one in multi_energy.cpp). We will have to go through
all the other code and change which one is being called.

***********************************************************/

/***********************************************************
Dependencies Begin Here!
***********************************************************/

/**********************************************************
 translates a BP (given as an integer) to two integerBases
**********************************************************/
function BP2_2(var bp, var base_i, var base_j)
{
   if (bp == 0) {base_i=0; base_j=3;}
   else if (bp == 1) {base_i=1; base_j=2;}
   else if (bp == 2) {base_i=2; base_j=1;}
   else if (bp == 3) {base_i=3; base_j=0;}
   else if (bp == 4) {base_i=2; base_j=3;}
   else if (bp == 5) {base_i=3; base_j=2;}
   else
   {
	  throw("Wronf BasepairNr.!");
	  return false;
   }
}

/**********************************************************
 translates an integerBase to a characterBase (nucleotide)
**********************************************************/
function int2char(base) {
	if(base == 0) return 'A';
	if(base == 1) return 'C';
	if(base == 2) return 'G';
	if(base == 3) return 'U';
	throw("There's something wrong in your int-base!");
}


/**********************************************************
 translates the two bases (integer) of a BP to an integer
 value for the pair (0,1,2,3,4,5)
**********************************************************/
function BP2int(b1, b2) {
	var kind_of_BP = 0;
	if(b1 == 0 && b2 == 3) kind_of_BP = 0;
	else if (b1 == 1 && b2 == 2) kind_of_BP = 1;
	else if (b1 == 2 && b2 == 1) kind_of_BP = 2;
	else if (b1 == 3 && b2 == 0) kind_of_BP = 3;
	else if (b1 == 2 && b2 == 3) kind_of_BP = 4;
	else if (b1 == 3 && b2 == 2) kind_of_BP = 5;
	else
	{
		throw("Wrong basepair (" + int2char(b1) + " - " + int2char(b2) + ") !");
	}
	return kind_of_BP;
}

/***********************************************************
 tests, whether one of three values is higher than MAX_DOUBLE,
 if so, the sum is set to MAX_DOUBLE as well
***********************************************************/

function Sum_MaxDouble3(var sum1, var sum2, var sum3)
{
   if ((sum1 == Number.MAX_VALUE) || (sum2 == Number.MAX_VALUE) || (sum3 == Number.MAX_VALUE))
	  return Number.MAX_VALUE;
   else
	  return sum1+sum2+sum3;
}



/***********************************************************
Dependencies End Here!
***********************************************************/

/**********************************************************************************
Global Constants                                                                   
**********************************************************************************/

var RT = 0.616;
var terminalAU = 0.5;
var loop_destabilizing_energies = new Array[
/***********************************************************
	 INTERNAL     BULGE     HAIRPIN     SIZE
  -------------------------------------------------------
***********************************************************/
	  0.0,      3.80,       0.0,     //  1  
	  0.0,      2.80,       0.0,     //  2
	  0.0,      3.20,      5.70,     //  3
	 1.70,      3.60,      5.60,     //  4
	 1.80,      4.00,      5.60,     //  5
	 2.00,      4.40,      5.40,     //  6
	 2.20,      4.59,      5.90,     //  7
	 2.30,      4.70,      5.60,     //  8
	 2.40,      4.80,      6.40,     //  9
	 2.50,      4.90,      6.50,     // 10
	 2.60,      5.00,      6.60,     // 11
	 2.70,      5.10,      6.70,     // 12
	 2.78,      5.19,      6.78,     // 13
	 2.86,      5.27,      6.86,     // 14 
	 2.94,      5.34,      6.94,     // 15
	 3.01,      5.41,      7.01,     // 16
	 3.07,      5.48,      7.07,     // 17
	 3.13,      5.54,      7.13,     // 18
	 3.19,      5.60,      7.19,     // 19
	 3.25,      5.65,      7.25,     // 20
	 3.30,      5.71,      7.30,     // 21
	 3.35,      5.76,      7.35,     // 22
	 3.40,      5.80,      7.40,     // 23
	 3.45,      5.85,      7.44,     // 24
	 3.49,      5.89,      7.49,     // 25
	 3.53,      5.94,      7.53,     // 26
	 3.57,      5.98,      7.57,     // 27
	 3.61,      6.02,      7.61,     // 28
	 3.65,      6.05,      7.65,     // 29
	 3.69,      6.09,      7.69      // 30
]

var stacking_energies = new Array[
//           Y                           Y                           Y                           Y 
//------------------------    ------------------------    ------------------------    ------------------------   
// A      C      G      U      A      C      G      U      A      C      G      U      A      C      G      U  
//------------------------    ------------------------    ------------------------    ------------------------   
//        5' --> 3'                   5' --> 3'                   5' --> 3'                   5' --> 3'     
//           AX                          AX                          AX                          AX 
//           AY                          CY                          GY                          UY 
//        3' <-- 5'                   3' <-- 5'                   3' <-- 5'                   3' <-- 5'    
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -0.90, 
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.20,   0.0,   
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.10,   0.0, -0.60, 
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.10,   0.0, -1.40,   0.0,    

//           Y                           Y                           Y                           Y 
//------------------------    ------------------------    ------------------------    ------------------------   
// A      C      G      U      A      C      G      U      A      C      G      U      A      C      G      U  
//------------------------    ------------------------    ------------------------    ------------------------   
//        5' --> 3'                   5' --> 3'                   5' --> 3'                   5' --> 3'     
//           CX                          CX                          CX                          CX 
//           AY                          CY                          GY                          UY 
//        3' <-- 5'                   3' <-- 5'                   3' <-- 5'                   3' <-- 5'    
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.10,   0.0,   0.0,   0.0,   0.0,   
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -3.30,   0.0,   0.0,   0.0,   0.0,   0.0,    
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.40,   0.0, -1.40,   0.0,   0.0,   0.0,   0.0,   
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.10,   0.0, -2.10,   0.0,   0.0,   0.0,   0.0,   0.0,    

//           Y                           Y                           Y                           Y 
//------------------------    ------------------------    ------------------------    ------------------------   
// A      C      G      U      A      C      G      U      A      C      G      U      A      C      G      U  
//------------------------    ------------------------    ------------------------    ------------------------   
//        5' --> 3'                   5' --> 3'                   5' --> 3'                   5' --> 3'     
//           GX                          GX                          GX                          GX 
//           AY                          CY                          GY                          UY 
//        3' <-- 5'                   3' <-- 5'                   3' <-- 5'                   3' <-- 5'    
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.40,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.30, 
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -3.40,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -2.50,   0.0,   
  0.0,   0.0,   0.0,   0.0,   0.0, -3.30,   0.0, -1.50,   0.0,   0.0,   0.0,   0.0,   0.0, -2.10,   0.0, -0.50,  
  0.0,   0.0,   0.0,   0.0, -2.20,   0.0, -2.50,   0.0,   0.0,   0.0,   0.0,   0.0, -1.40,   0.0,  1.30,   0.0,     

//           Y                           Y                           Y                           Y 
//------------------------    ------------------------    ------------------------    ------------------------   
// A      C      G      U      A      C      G      U      A      C      G      U      A      C      G      U  
//------------------------    ------------------------    ------------------------    ------------------------   
//        5' --> 3'                   5' --> 3'                   5' --> 3'                   5' --> 3'     
//           UX                          UX                          UX                          UX 
//           AY                          CY                          GY                          UY 
//        3' <-- 5'                   3' <-- 5'                   3' <-- 5'                   3' <-- 5'    
  0.0,   0.0,   0.0, -1.30,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.00,   0.0,   0.0,   0.0,   0.0,   
  0.0,   0.0, -2.40,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.50,   0.0,   0.0,   0.0,   0.0,   0.0,    
  0.0, -2.10,   0.0, -1.00,   0.0,   0.0,   0.0,   0.0,   0.0, -1.40,   0.0,  0.30,   0.0,   0.0,   0.0,   0.0,   
-0.90,   0.0, -1.30,   0.0,   0.0,   0.0,   0.0,   0.0, -0.60,   0.0, -0.50,   0.0,   0.0,   0.0,   0.0,   0.0
]

/**********************************************************
 bulge_energy.cpp
**********************************************************/


function BulgeEnergy(size, bp_i, bp_j, bp_before) {
	
	var energy = 0.0;
	var base_i_before, base_j_before;

	BP2_2(bp_before, base_i_before, base_j_before);

	/**********************************************************
	 No bulges of size 0:
	**********************************************************/
	if(size == 0)
	{
		throw("No bulge!");
	}

	/**********************************************************
	 Loop destabilizing energies
	**********************************************************/
	if(size <= 30)
		energy += loop_destabilizing_energies[3*size-2];
	else
	{
		energy += loop_destabilizing_energies[3*30-2];
		energy += 1.75*RT*Math.log(size/30.0);
	}

	/**********************************************************
	 Additional bulge-energy (size 1 : stacking, >1 : nothing but term, A-U penalty):
	**********************************************************/
	if (size == 1)
		energy += stacking_energies[64*bp_i+16*base_i_before+4*bp_j+base_j_before];
	else
	{
		if ((bp_before == 0) || (bp_before == 3) || (bp_before == 4 || bp_before == 5))
			energy += terminalAU;
		if ((BP2int(bp_i, bp_j) == 0) || (BP2int(bp_i, bp_j) == 3) || (BP2int(bp_i, bp_j) == 4) || (BP2int(bp_i, bp_j) == 5))
			energy += terminalAU;
	}
	return energy;
}

/**********************************************************
 stacking_energy.cpp
**********************************************************/


function StackingEnergy(bp_i, bp_j, bp_before) {

	var energy = 0.0;
	var base_I_before, base_j_before;

	BP2_2(bp_before, base_i_before, base_j_before);
	energy += stacking_energies[64*bp_i+16*base_i_before+4*bp_j+base_j_before];

	return energy;
}

/**********************************************************
 end_energy.cpp
 **********************************************************/

 /*******************************************************
 identifies connections concerning the closing base pairs
 of the stems and concerning the free bases adjacent to them

 stem_num = number of stems
 vorgaenger = contains all stem ends and their positions in BP_Order
 BaseConnections = contains row-wise the connections concerning bases
 BasePairConnections = contains row-wise the connections concerning base pairs
*******************************************************/
function EndConnections(stem_num, vorgaenger, EndBaseConnections, EndBasePairConnections) {
	var bases_between_stems = new Array(stem_num);	//contains the number of bases between the stems (in each coord. the bases "before" the
													//stem, only at the last coord. of bases_between_stems the free bases after
													//the last stem and before the closing BP are stored
	
	var one_connection = new Array(); //help vector for one connection
	var i;

	// number of free bases between the stem-ending base pairs
	for(i = 0; i < stem_num-1; i++)
		bases_between_stems[i] = BP_Order[vorgaenger[i+1]][0] - BP_Order[vorgaenger[i]][1] - 1;
	
	bases_between_stems[stem_num-1] = brackets.length - BP_Order[vorgaenger[stem_num-1]][1] - 1;

	// finding all connections: (concerning free bases)
	//**************************************************
	for(i = 0; i < stm_num-1; i++) {
		// if there is only one free base between the stems and if there is at least one free base after the next stem
		if(bases_between_stems[i] >= 1) {
			one_connection.push(BP_Order[vorgaenger[i]][1] + 1);
			if(bases_between_stems > 1) {
				EndBaseConnections.push(one_connection);
				one_connection.length = 0; // practically the same as .clear();
				one_connection.push(BP_Order[vorgaenger[i+1]][0] - 1);
			}
		} else {
			EndBaseConnections.push(one_connection);
			one_connection.length = 0; 
		}
	}

	if(bases_between_stems[stem_num - 1] >= 1) {
		one_connection.push(BP_Order[vorgaenger[stem_num-1]][1]+1);
		EndBaseConnections.push(one_connection);
		one_connection.length = 0; 
	} else {
		EndBaseConnections.push(one_connection);
		one_connection.length = 0;
	}

	// finding all connections: (concerning stems)
	//*************************************************
	one_connection.push(vorgaenger[0]);

	for(i = 0; i < stem_num-1; i++) {
		if(bases_between_stems[i] == 1)
			one_connection.push(vorgaenger[i+1]);
		else {
			EndBasePairConnections.push(one_connection);
			one_connection.length = 0;
			one_connection.push(vorgaenger[i+1]);
		}
	}
	EndBasePairConnections.push(one_connection);
	one_connection.length = 0;
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

function EndConnections(vorgaenger, EndBaseConnections, EndBaseSizes, EndBasePairConnections, EndBasePairSizes) {
	var stem_num = Math.max(BP_Order[numBP][3], 1);
	var num_bp_con = 0; 	//current number of BasePairConnections
	var num_base_con = 0;	//current number of BaseConnections
	var bp, base;	//current position in the current connection
	var bases_between_stems = new Array(stem_num);	//contains the number of bases between the stems (in each coord. the bases "before" the
													//stem, only at the last coord. of bases_between_stems the free bases after
													//the last stem and before the closing BP are stored

	var i;                                            

	// number of free bases between the stem-ending base pairs
	for(i = 0; i < stem_num-1; i++)
		bases_between_stems[i] = BP_Order[vorgaenger[i+1]][0] - BP_Order[vorgaenger[i]][1] - 1;
	bases_between_stems[stem_num-1] = struct_len - BP_Order[vorgaenger[stem_num-1]][1] - 1;

	// finding all connections: (concerning free bases)
	//**************************************************
	bp = 0;
	base = 0;

	for(i = 0; i < stem_num-1; i++) {

		// if there is only one free base between the stem and if there is atleast one free base after the next stem
		if(bases_between_stems[i] >= 1) {
			EndBaseConnections[num_base_con][base++] = BP_Order[vorgaenger[i]][1] + 1;
			if(bases_between_stems[i] > 1) {
				EndBaseSizes[num_base_con] = base;
				num_base_con++;
				base = 0;
				EndBaseConnections[num_base_con][base++] = BP_Order[vorgaenger[i+1]][0] - 1;
			}
		} else {
			EndBaseSizes[num_base_con] = base;
			num_base_con++;
			base = 0;
		}
	}

	if(bases_between_stems[stem_num - 1] >= 1) {
		EndBaseConnections[num_base_con][base++] = BP_Order[vorgaenger[stem_num-1]][1] + 1;
		EndBaseSizes[num_base_con] = base;
		num_base_con++;
		base = 0;
	} else {
		EndBaseSizes[num_base_con] = base;
		num_base_con++;
		base = 0;
	}

	// finding all connections: (concerning stems)
	//*************************************************
	EndBasePairConnections[num_bp_con][bp++] = vorgaenger[0];

	for(i = 0; i < stem_num-1; i++) {
		if(bases_between_stems[i] == 1)
			EndBasePairConnections[num_bp_con][bp++] = vorgaenger[i+1];
		else {
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

function BestEndConnectionEnergy(base_connection, pair_connection, bp_at_pair_connection) {
	var energy = 0.0;
	var min = Number.MAX_VALUE;
	var base_size = base_connection.length;
	var pair_size = pair_connection.length;
	var base_energy;
	var i, j;

	for(i = 0; i < base_size; i++) {
		min = Number.MAX_VALUE;

		/*There are 4 cases:  (a) BP b BP b
							  (b) b BP b BP b
							  (c) b BP b BP
							  (d) BP b BP*/

		if(i - 1 >= 0) { /* case BP_ORDER[pair_connection[i-1]][0] not possible */
			// usual case : base is located right of (= after) the stem
			// corresponds to cases (b) and (c) (from the second free base on)
			if(BP_ORDER[pair_connection[i-1]][1]+1 == base_connection[i])
				for(j = 0; j < 4; j++) {
					// energy and penalty
					base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j], BasePenalty(base_connection[i], j));
				}
		}

		if(i < pair_size) {
			// usual case : base is located left of (= before) the stem
			// corresponds to cases (b) and (c) (free bases in front of a following stem)
			if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
				for(j = 0; j < 4; j++) {
					base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][1] + 4*bp_at_pair_connection[i][0]+j], BasePenalty(base_connection[i],j));
					if(base_energy < min)
						min = base_energy;
				}

			// usual case : base is located right of (=after) the stem
			// corresponds to cases (a) and (d) (free bases after the stem)
			if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
				for(j = 0; j < 4; j++) {
					base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j], BasePenalty(base_connection[i], j));
					if (base_energy < min)
						min = base_energy;
				}
		}

		if (i + 1 < pair_size) {
			// usual case : base is located left of (=before) the stem
			// corresponds to cases (a) and (d) (free base in front of the stem)
			if(BP_Order[pair_connection[i+1]][0] -1 == base_connection[i])
				for(j = 0; j < 4; j++) {
					base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][1]+4*bp_at_pair_connection[i+1][0]+j], BasePenalty(base_connection[i], j));
					if(base_energy < min)
						min = base_energy;
				}
		}
		energy = Sum_MaxDouble(energy, min);
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

function EndConnectionEnergy(base_connection, base_size, pair_connection, pair_size, int_seq) {

	var energy = 0.0;
	var min;
	var i;

	for(i = 0; i < base_size; i++) {
		min = Number.MAX_VALUE;

		/*There are 4 cases:  (a) BP b BP b
							  (b) b BP b BP b
							  (c) b BP b BP
							  (d) BP b BP*/

		if(i - 1 >= 0) {/* case BP_Order[pair_connection[i-1]][0] not possible */
			// usual case : base is located right of (= after) the stem
			// corresponds to cases (b) and (c) (from the second free base on)
			if(BP_Order[pair_connection[i-1]][1] + 1 == base_connection[i])
				if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i-1]][1]]+4*int_seq[BP_Order[pair_connection[i-1]][0]]+int_seq[base_connection[i]]] < min)
					min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i-1]][1]]+4*int_seq[BP_Order[pair_connection[i-1]][0]]+int_seq[base_connection[i]]];
		}

		if( i < pair_size) {
			// usual case : base is located left of (= before) the stem
			// corresponds to cases (b) and (c) (free bases in front of a following stem)
			if(BP_Order[pair_connection[i]][0]-1 == base_connection[i])
				if (single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]] < min)
					min = single_base_stacking_energy[64+16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]];

			// usual case : base is located right of (= after) the stem
			// corresponds to cases (a) and (d) (free bases after the stem)
			if(BP_Order[pair_connection[i]][1]+1 == base_connection[i])
				if (single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]] < min)
					min = single_base_stacking_energy[16*int_seq[BP_Order[pair_connection[i]][1]]+4*int_seq[BP_Order[pair_connection[i]][0]]+int_seq[base_connection[i]]];

		}

		if(i+1 < pair_size) {
			// usual case : base is located left of (= before) the stem
			// corresponds to cases (a) and (d) (free base in front of the stem)
			if(BP_Order[pair_connection][i+1]][0]-1 == base_connection[i])
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

function externBestEnergy() {
	// in the external loop we have to take into account the ending stems, the free base between them, and the dangling ends
	/*e.g.:    bp_i - bp_j                    bp_i - bp_j
		  closing_i - closing_j freebase closing_i - closing_j free_base*/

	//finding all predecessors of the external loop:
	/*checking the structure vector (brackets), start at the first opening bracket, jump to its closing counterpart (or to the following pos.),
	  number of stems in the external loop++, then looking for the next opening bracket, jump to its closing counterpart +1, ....*/

	var bp_i, bp_j;
	var i, a, j, c;
	var vorg_last = new Array(BP_Order[numBp][3]);
	for(a = 0; a < BP_Order[numBp][3]; a++)
		vorg_last[a] = -1;

	var stems = 0;
	var pos = 0;
	while(pos < struct_len) {
		if(brackets[pos] == '(') {
			vorg_last[stems] = BP_Pos_Nr[pos];
			stems++;
			pos = BP_Order[BP_Pos_Nr[pos]][1] + 1;
		} else {
			pos++;
		}
	}

	// initially the last row of D is filled with values of the row before, if there are no dang-ends, no further changes have to be done to the last row
	var bp_assign, vg, free_base;
	for(bp_assign = 0; bp_assign < 6; bp_assign++)
		D[numBP][bp_assign] = D[numBP-1][bp_assign];

	var min_dang;
	var front_dang = 0.0;
	var base_energy;

	var EndBaseConnections; // all connections, connecting free bases
	var EndBasePairConnections; // all connections, concerning stems (BP)
	EndConnections(Math.max(BP_Order[numBP][3], 1), vorg_last, EndBaseConnections, EndBasePairConnections);

	for(bp_assign = 0; bp_assign < 6; bp_assign++) {
		for(vg = 0; vg < BP_Order[numBP][3]; vg++)
			Trace[numBP][bp_assign][vg][0] = vorg_last[vg];
		// during the traceback go to the same assignment as this position in the artificial last row of D and Trace
		Trace[numBP][bp_assign][0][1] = bp_assign; //--> this can be overwritten later, but if the second if-loop is not reached, this iniitialization is needed

		// if there are free bases in front of the first opening bracket
		if (BP_Order[numBP-1][0] > 0) {
			min_dang = Number.MAX_VALUE;
			BP2_2(bp_assign, bp_i, bp_j);
			for(free_base = 0; free_base < 4; free_base++) {
				base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_j+4*bp_i+free_base], BasePenalty(BP_Order[numBP-1][0]-1,free_base));
				if(min_dang > base_energy)
					min_dang = base_energy;
			}
			front_dang = min_dang;
		}

		// if there is a structural fragment after the closing bracket of the last BP:
		if(BP_Order[numBP-1][1] < strct_len-1) {
			var min = Number.MAX_VALUE;
			var energy = 0.0;
			var energy_help = 0.0, energy_unbound, energy_stems;
			var pair_size, base_size;

			var bp_at_pair_connection; // base pair assignment for a connection
			var MIN_bp_at_pair_connection; // base pair assignment for a connection with min. free energy

			for(i = 0; i < EndBaseConnections.length; i++) { // EndBaseConnections.length = EndBasePairConnections.length
				pair_size = EndBasePairConnections[i].length;
				base_size = EndBaseConnections[i].length;

				bp_at_pair_connection = new Array(pair_size);
				MIN_bp_at_pair_connection = new Array(pair_size);
				for(a = 0; a < pair_size; a++) {
					bp_at_pair_connection[a] = new Array(2);
					MIN_bp_at_pair_connection[a] = new Array(2);
					for(b = 0; b < 2; b++) {
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

				if(EndBasePairConnections[i][0] == numBP-1) {
					start_pos = 1;
					BP2_2(bp_assign, bp_i, bp_j);
					bp_at_pair_connection[0][0] = bp_i;
					bp_at_pair_connection[0][1] = bp_j;
					start_power = pair_size - 2;
				} else {
					start_pos = 0;
					start_power = pair_size - 1;
				}

				min = Number.MAX_VALUE;
				var int_anz;
				for(int_anz = 0; int_anz < Math.pow(6, start_power+1); int_anz++) {
					var power = start_power;
					var value = int_anz;
					var pos = start_pos;
					var bp_help;

					// create an assignment
					while(power >= 0) {
						bp_help = value / Math.pow(6, power);

						BP2_2(bp_help, bp_i, bp_j);
						bp_at_pair_connection[pos][0] = bp_i;
						bp_at_pair_connection[pos][1] = bp_j;

						value = fmod(value, Math.pow(6, power));
						power--;
						pos++;
					}

					// energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here, energy of the stem will be added later)
					energy_unbound = BestEndConnectionEnergy(EndBaseConnections[i], EndBasePairConnections[i], bp_at_pair_connection);

					// add the energy of the stems (depending on the closing BP), BUT: not closing BP!
					energy_stems = 0.0;
					for(j = 0; j < pair_size; j++)
						energy_stems = Sum_MaxDouble(energy_stems, D[EndBasePairConnections[i][j]][BP2int(bp_at_pair_connection[j][0], bp_at_pair_connection[j][1])]);

					energy_help = Sum_MaxDouble(energy_stems, energy_unbound);

					if(min > energy_help) {
						min = energy_help;
						for(a = 0; a < pair_size; a++) {
							MIN_bp_at_pair_connection[a][0] = bp_at_pair_connection[a][0];
							MIN_bp_at_pair_connection[a][1] = bp_at_pair_connection[a][1];
						}
					}
				}// for int_anz

				// if no valid assignment was found, set a random one (this is the case, if in the external loop a BP exists that has no valid assignment
				// concerning the constraints --> usually this error is filtered before)
				if(min == Number.MAX_VALUE)
					for(a = 0; a < pair_size; a++) {
						var rand_bp = RandomBasePair();
						BP2_2(rand_bp, MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);
					}

				energy = Sum_MaxDouble(energy, min);

				// store best assignment of the predecessors in the traceback
				for(a = 0; a < pair_size; a++)
					for(vg = 0; vg < BP_Order[numBP][3]; vg++)
						if(EndBasePairConnections[i][a] == Trace[numBP][bp_assign][vg][0])
							Trace[numBP][bp_assign][vg][1] = BP2int(MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);

			} // for i

			D[numBP][bp_assign] = Sum_MaxDouble(energy, front_dang);
		} // if

	} // for bp_assign

	for(i = 0; i < EndBaseConnections.length; i++)
		EndBaseConnections[i].length = 0;
	EndBaseConnections.length = 0;

	for(i = 0; i < EndBasePairConnections.length; i++)
		EndBasePairConnections[i].length = 0;
	EndBasePairConnections.length = 0;
}

/****************************************************************************************

 finds energy of the external loop (for a given sequence!!)

****************************************************************************************/

function externEnergy(int_seq) {

   // in the external loop we have to take into account the ending stems, the free base between them, and the dangling ends
   /*e.g.:    bp_i - bp_j                    bp_i - bp_j
		 closing_i - closing_j freebase closing_i - closing_j free_base*/

  var bp_i, bp_j, a, free_base;
  var i, j;

  var stems = 0;
  var pos = 0;

  var EndBaseConnections; // all connections, concerning free bases
  var EndBasePairConnections; // all connections, concerning stems
  var EndBaseSizes; // sizes of the base connections
  var EndBasePairSizes; // sizes of the stem connections
  var vorg_last; // BP_Order-pos. of the BPs that are predecssors of the external loop

  var num_conns;
  var energy = 0.0;
  var energy_unbound;
  var pair_size, base_size;

  var con_size = BP_Order[numBP][3] + 2; // in order that base and basepairconn. have the same size, no conn. can be bigger than #stems+2 and all together
										 // there can not be more (since only adjacent free base are taken into account)

   //finding all predecessors of the external loop:
   /*checking the structure vector (brackets), start at the first opening bracket, jump to its closing counterpart (or to the following pos.),
	 number of stems in the external loop++, then looking for the next opening bracket, jump to its closing counterpart +1, ....*/
  vorg_last = new Array(BP_Order[numBP][3]);
  for(a = 0; a < BP_Order[numBP][3]; a++)
	vorg_last[a] = -1;

  while(pos < struct_len) {
	if(brackets[pos] == '(') {
	  vorg_last[stems++] = BP_Pos_Nr[pos];
	  pos = BP_Order[BP_Pos_Nr[pos]][1] + 1;
	} else {
	  pos++;
	}
  }

  // if there are free bases in front of the first opening brackets
  if(BP_Order[numBP-1][0] > 0) {
	bp_i = int_seq[BP_Order[numBP-1][0]];
	bp_j = int_seq[BP_Order[numBP-1][1]];
	free_base = int_seq[BP_Order[numBP-1][0]-1];
	energy = Sum_MaxDouble3(energy, single_base_stacking_energy[64+16*bp_j+4*bp_i+free_base], BasePenalty(BP_Order[numBP-1][0]-1, free_base));
	// if there are more free bases in front of the first opening bracket, the penalties are considered seperatly in get_basePart_Energyi n
	// multi_energy.cpp
  }                    

  // since these arrays here cannot be created dynamically, an array of maximal size is created
  // this is not true in javascript, but i left it like this
  EndBaseConnections = new Array(con_size);
  EndBaseSizes = new Array(con_size);
  for(i = 0; i < con_size; i++) {
	EndBaseConnections[i] = new Array(con_size);
	EndBaseSizes[i] = -1;
	for(j = 0; j < con_size; j++)
	  EndBaseConnections[i][j] = -1;
  }

  EndBasePairConnections = new Array(con_size);
  EndBasePairSizes = new Array(con_size);
  for(i = 0; i < con_size; i++) {
	EndBasePairConnections[i] = new Array(con_size);
	EndBasePairSizes[i] = -1;
	for(j = 0; j < con_size; j++)
	  EndBasePairConnections[i][j] = -1;
  }

  EndConnections(vorg_last, EndBaseConnections, EndBaseSizes, EndBasePairConnections, EndBasePairSizes);

  // count number of connections
  num_conns = 0;
  for(i = 0; i < Math.Max(BP_Order[numBP][2], BP_Order[numBP][3] + 2); i++) {
	if( (EndBaseConnections[i][0] != -1) || (EndBasePairConnections[i][0] != -1))
	  num_conns++;
	else
	  break;
  }

  for(i = 0; i < num_conns; i++) {
	pair_size = EndBasePairSizes[i];
	base_size = EndBaseSizes[i];

	// energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here, energy of
	// the stem will be added later)
	energy_unbound = EndConnectionEnergy(EndBaseConnections[i], base_size, EndBasePairConnections[i], pair_size, int_seq);
	energy = Sum_MaxDouble(energy, energy_unbound);
  }

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

function EndConnectionBestFreeBases(base_connection, pair_connection, bp_at_pair_connection) {
  var min = Number.MAX_VALUE;
  var base_energy;
  var base_size = base_connection.length;
  var pair_size = pair_connection.length;
  var min_base = 0;
  var min_bases = new Array(base_size);
  var i, j;
  for(i = 0; i < base_size; i++)
	min_bases[i] = -1;

  for(i = 0; i < base_size; i++) {
	//for each base: testing whether it is adjacent to a stem or two, minimize the energy
	min = Number.MAX_VALUE;

	if(i - 1 >= 0) /* case BP_Order[pair_connection[i-1]][0] not possible */ {
	  // usual case : base is located right of (= after) the stem
	  if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
		for(j = 0; j < 4; j++) {
		  base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j], BasePenalty(base_connection[i], j));
		  if(base_energy < min) {
			min = base_energy;
			min_base = j;
		  }
		}
	}

	if(i < pair_size) {
	  // usual case : base is located left of (= before) the stem
	  if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
		for(j = 0; j < 4; j++) {
		  base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
		  if(base_energy < min) {
			min = base_energy;
			min_base = j;
		  }
		}

	  // usual case : base is located right of (= after) the stem
	  if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
		for(j = 0; j < 4; j++) {
		  base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
		  if(base_energy < min) {
			min = base_energy;
			min_base = j;
		  }
		}
	}

	if (i+1 < pair_connection.length) {
	  // usual case : base is located left of (= before) the stem
	  if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
		for(j = 0; j < 4; j++) {
		  base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][1]+4*bp_at_pair_connection[i+1][0]+j],BasePenalty(base_connection[i],j));
		  if(base_energy < min) {
			min = base_energy;
			min_base = j;
		  }
		}
	}

	min_bases[i] = min_base;

  }

  return min_bases;
}



/**********************************************************
 hairpin_energy.cpp
 **********************************************************/

 /******************************************************
identifies the energy boni for special tetra loops
******************************************************/

function tetra_loop_energy(str) {
  if (str === "GGGGAC") return -3.00;
  if (str === "GGUGAC") return -3.00;
  if (str === "CGAAAG") return -3.00;
  if (str === "GGAGAC") return -3.00;
  if (str === "CGCAAG") return -3.00;
  if (str === "GGAAAC") return -3.00;
  if (str === "CGGAAG") return -3.00;
  if (str === "CUUCGG") return -3.00;
  if (str === "CGUGAG") return -3.00;
  if (str === "CGAAGG") return -2.50;
  if (str === "CUACGG") return -2.50;
  if (str === "GGCAAC") return -2.50;
  if (str === "CGCGAG") return -2.50;
  if (str === "UGAGAG") return -2.50;
  if (str === "CGAGAG") return -2.00;
  if (str === "AGAAAU") return -2.00;
  if (str === "CGUAAG") return -2.00;
  if (str === "CUAACG") return -2.00;
  if (str === "UGAAAG") return -2.00;
  if (str === "GGAAGC") return -1.50;
  if (str === "GGGAAC") return -1.50;
  if (str === "UGAAAA") return -1.50;
  if (str === "AGCAAU") return -1.50;
  if (str === "AGUAAU") return -1.50;
  if (str === "CGGGAG") return -1.50;
  if (str === "AGUGAU") return -1.50;
  if (str === "GGCGAC") return -1.50;
  if (str === "GGGAGC") return -1.50;
  if (str === "GUGAAC") return -1.50;
  if (str === "UGGAAA") return -1.50;
  return 0.0;
}

/******************************************************************

finds best boni for tetraloops depending on the closing BP + the
first mismatch
	 ==> fct. is not longer used when dealing with seq. constraints
*******************************************************************/
function best_tetra_loop_energy(bp_i, bp_j, i, j) {

  if ((bp_i == 0) && (bp_j == 3)) {
	if ((i==2) && (j==0)) return = -2.00;
  }
  else if ((bp_i == 1) && (bp_j == 2)) {
	if ((i==2) && (j==0)) return -3.00;
	else if ((i==2) && (j==2)) return -2.5;
	else if ((i==3) && (j==1)) return -2.0;
	else if ((i==3) && (j==2)) return -3.0;
  }
  else if ((bp_i == 2) && (bp_j == 1)) {
	if ((i==2) && (j==0)) return -3.00;
	else if ((i==2) && (j==2)) return -1.5;
	else if ((i==3) && (j==0)) return -1.5;
  }
  else if ((bp_i == 3) && (bp_j == 0)) {
	if ((i==2) && (j==0)) return -1.50;
  }
  else if ((bp_i == 3) && (bp_j == 2)) {
	if ((i==2) && (j==0)) return -2.50;
  }

  return 0.0;
}

/*************************************************************

calculation of the minimal free energy for an HL with given
size and given closing BP

*************************************************************/
function BestHairpinLoopEnergy(bp_pos, size, bp_i, bp_j) {
  var energy = 0.0;
  var min = Number.MAX_VALUE;
  var energy_help;

  var bp_pos_i = BP_Order[bp_pos][0];
  var bp_pos_j = BP_Order[bp_pos][1];

  var i, j, i2, j2;

  // loop destabilizing energies
  if( size <= 30)
	energy += loop_destabilizing_energies[3*size-1];
  else {
	energy += loop_destabilizing_energies[3*30-1];
	energy += 1.75*RT*Math.log(size/30.0);
  }

  // terminal mismatch (size > 3)
  if (size > 3) {
	  // The exact assignment of the interior loop bases (not adjacent to the closing BP) is not important here.
	  // Their valid assignments can be restricted by the constraints, but one assignment per position has to be
	  // valid und this one is chosen for the best energy.

	  // If the HPLoop would only consist of C's, a penalty of 0.3*size+1.6 has to be added.
	  // This means that just in case that all bases in the loop are restricted to C by the
	  // constraints the energy has to be corrected.
	  // NOT TAKEN INTO ACCOUNT, SINCE THIS IS NOT TAKEN INTO ACCOUNT IN THE VIENNA PACKAGE AS WELL.

	  // bool onlyCs = true;
	  // for (int p=bp_pos_i+1; p<bp_pos_j; p++)
	  //    if ( !((seq_constraints[p][1] == 1) && (Sum_SeqConst(p) == 1)))
	  //    {
	  //       onlyCs = false;
	  //       break;
	  //    }

	  // if (onlyCs)
	  //    energy = Sum_MaxDouble(energy,0.3*size+1.6);

	  if(size != 4) {
		// finding the minimal terminal mismatch energy for the given closing BP
		// the remaining bases in the loop don't contribtue to the energy, i.e. they
		// are not fixed until the traceback, their penalty is also considered there
		// (since at least one base has to be valid for each position, the exact
		// assignment has not to be fixed here)
		min = Number.MAX_VALUE;
		for(i= 0; i < 4; i++) {
		  for(j = 0; j < 4; j++) {
			energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1, i), BasePenalty(bp_pos_j-1, j), mismatch_energies_hairpin[64*bp_i+16*i+4*bp_j+j]);
			if(energy_help < min)
			  min = energy_help;
		  }
		}
		energy = Sum_MaxDouble(energy, min);
	  } else {
		var tetra_plus_closing = new Array(6);
		tetra_plus_closing[0] = int2char(bp_i);
		tetra_plus_closing[5] = int2char(bp_j);

		 // consider all cases of the tetraloop (4*4*4*4) and add term-mismatch and
		 // possibly a bonus, furthermore add the penalties for all 4 bases
		 //      i2   j2
		 //   i           j
		 //     bp_i-bp_j
		 //       .   .

		 min = Number.MAX_VALUE;
		 for(i = 0; i < 4; i++)
		  for(j = 0; j < 4; j++)
			for(i2 = 0; i2 < 4; i2++)
			  for(j2 = 0; j2 < 4; j2++) {
				energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1, i), BasePenalty(bp_pos_i+2, i2), BasePenalty(bp_pos_j-2, j2), BasePenalty(bp_pos_j-1, j));
				tetra_plus_closing[1] = int2char(i);
				tetra_plus_closing[2] = int2char(i2);
				tetra_plus_closing[3] = int2char(j2);
				tetra_plus_closing[4] = int2char(j);
				energy_help = Sum_MaxDouble3(energy_help, tetra_loop_energy(tetra_plus_closing), mismatch_energies_hairpin[64*bp_i+16*i+4*bp_j+i]);

				if(energy_help < min)
				  min = energy_help;

			  }
		 energy = Sum_MaxDouble(energy, min);        
	  }

  } else if(size == 3) {
	min = Number.MAX_VALUE;
	for(i = 0; i < 4; i++)
	  for(j = 0; j < 4; j++)
		for(i2 = 0; i2 < 4; i2++) {
		  energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1, i), BasePenalty(bp_pos_i+2, i2), BasePenalty(bp_pos_j-1, j));

		  // If the HL would only consist of Cs, a penalty of 1.4 has to be added.
		  // This is just the best solution, if all other bases are forbidden by the constraints.
		  if((i == 1) && (i2==1) && (j==1))
			energy_help = Sum_MaxDouble(energy_help, Ctriloop);

		  // If the HL would only consist of Gs, a term of -2.2 has to be added
		  // This is not considered in the ViennaPackage, so do we
		  // if ((i == 2) && (i2 == 2) && (j ==2))
		  //  energy_help = Sum_MaxDouble(energy_help, Gtriloop);

		  if(energy_help < min)
			min = energy_help;

		}

  energy = Sum_MaxDouble(energy, min);

  if(((bp_i == 0) && (bp_j == 3)) || ((bp_i == 3) && (bp_j == 0)) || ((bp_i == 2) && (bp_i==3)) || ((bp_i==3) && (bp_j==2)))
	energy = Sum_MaxDouble(energy, terminalAU);
	
  }
  else
  {
	throw "HairpinLoopSize too small!";
  }
  return energy;
}

/*************************************************************

calculates energy fractions of a HL that are important for 
the difference of two assignments of it
(size and closing BP are given)

*************************************************************/
function HairpinLoopEnergy(size, pos_i, int_seq) {
  var energy = 0.0;
  var i, k;
  var tetra_loop_assign;
  var pos_j = pos_i + size + 1;

  var bp_i = int_seq[pos_i];
  var bp_j = int_seq[pos_j];
  var mis_i = int_seq[pos_i+1];
  var mis_j = int_seq[pos_j - 1];

  var energy_help = 0.0;

  // loop_destabilizing_energies (not necessary here, since this function is
  // just used for the difference of two assignments of the same HL)
  //------------------------------
  /*if (size <= 30)
	energy += loop_destabilizing_energies[3*size-1];
  else
  {
	energy += loop_destabilizing_energies[3*30-1];
	energy += 1.75*RT*log((double)(size/30.0));
  }*/

  // terminal mismatch (for size > 3)
  if (size > 3) {

	  //**************************************************************
	  // If the HPLoop would only consist of C's, a penalty of 0.3*size+1.6 has to be added.
	  // NOT TAKEN INTO ACCOUNT, SINCE THIS IS NOT TAKEN INTO ACCOUNT IN THE VIENNA PACKAGE AS WELL.
	  /*only_Cs = 1;
	  for (i=1; i<=size; i++)
		 if (int_seq[pos_i+i] != 1)
		 {
			only_Cs = 0;
			break;
		 }
	  if ((only_Cs == 1) && (size > 3))
		 energy += 0.3*size+1.6;*/
	  //**************************************************************

	  // penalties for all bases of the loop
	  energy_help = 0.0;
	  for(i = 1; i <= size; i++)
		energy_help = Sum_MaxDouble(energy_help, BasePenalty(pos_i+i, int_seq[pos_i+1]));

	  // terminal mismatch energy for the closing BP and its penalty
	  energy = Sum_MaxDouble4(energy, energy_help, mismatch_energies_hairpin[64*bp_i+16*mis_i+4*bp_j+mis_j], PairPenalty(pos_i, pos_j, bp_i, bp_j));

	  if(size == 4) {
		// for tetraloops consider special term
		tetra_loop_assign = new Array(6);
		for(k = 0; k < 6; k++)
		  tetra_loop_assign[k] = int2char(int_seq[pos_i+k]);
		energy = Sum_MaxDouble(energy, tetra_loop_energy(tetra_loop_assign));
	  }
  } else if(Size == 3) {
	energy_help = Sum_MaxDouble3(BasePenalty(pos_i+1, mis_i), BasePenalty(pos_i+2, int_seq[pos_i+2]), BasePenalty(pos_j-1, mis_j));
	energy = Sum_MaxDouble3(energy, energy_help, PairPenalty(pos_i, pos_j, bp_i, bp_j));
  
	if(((bp_i==0) && (bp_j == 3)) || ((bp_i ==3) && (bp_j == 0)) || ((bp_i == 2) && (bp_j==3)) || ((bp_i==3) && (bp_j==2)))
	  energy = Sum_MaxDouble(energy, terminalAU);

	// If the HL would only consist of Cs, a penalty of 1.4 has to be added
	if((int_seq[pos_i+1] == 1) && (int_seq[pos_i+2] == 1) && (int_seq[pos_i+3] == 1))
	  energy = Sum_MaxDouble(energy, Ctriloop);

	// If the HL would only consist of Gs a term of -2.2 has to be added
	// This is not considered in the ViennaPackage so do we
	// if ((int_seq[pos_i+1] == 2) && (int_seq[pos_i + 2] == 2) && (int_seq[pos_i+3] == 2))
	//    energy = Sum_MaxDouble(energy, Gtriloop);

  } else {
	throw("HairpinLoopSize too small");
  }

  // if closingHL the last BP in a stack
  energy = Sum_MaxDouble(energy, Zero_or_StemEndAU(BP_Pos_Nr[pos_i], BP2int(bp_i, bp_j)));

  return energy;
}

/************************************************************

interior_energy.cpp

************************************************************/


/************************************************************

finds the minimal energy of an interior loop, if the
assignments of both closing BPs are known

************************************************************/

function BestInteriorLoopEnergy(bp_pos, leftSize, rightSize, bp_i, bp_j, bp_before) {
   var energy = 0.0;
   var asym = 0.0;
   var min = Number.MAX_VALUE;
   var energy_help;

   var size = leftSize + rightSize;
   var bp_new = BP2int(bp_i, bp_j);
   var base_i_before, base_j_before;
   var bp_pos_i, bp_pos_j;

   BP2_2(bp_before,base_i_before,base_j_before);
   bp_pos_i = BP_Order[bp_pos][0];
   bp_pos_j = BP_Order[bp_pos][1];

   var x, y, z, x1, x2, y1, y2, i1, i2, i3, i4;

   // no bulges and loops of size 0:
   //******************************************
   if ((leftSize == 0) || (rightSize == 0)) {
	  throw("No interior loop!");
   }

   // special cases:
   //******************************************
   if ((leftSize == 1) && (rightSize == 1)) {
	  for (x=0; x<4; x++)
		 for (y=0; y<4; y++) {
			energy_help = Sum_MaxDouble3(interior_loop_1_1_energy[96*bp_new+24*x+4*bp_before+y], BasePenalty(bp_pos_i+1,x), BasePenalty(bp_pos_j-1,y));
			if (energy_help < min)
			   min = energy_help;
		 }
	  energy = min;
   } else if ((leftSize == 1) && (rightSize == 2)) {
	  for (x=0; x<4; x++)
		 for (y=0; y<4; y++)
			for (z=0; z<4; z++) {
			   energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,x), BasePenalty(bp_pos_j-1,y), BasePenalty(bp_pos_j-2,z));
			   energy_help = Sum_MaxDouble(energy_help, interior_loop_1_2_energy[384*bp_new+96*z+24*x+4*bp_before+y]);
			   if (energy_help < min)
				  min = energy_help;
			}
	  energy = min;
   } else if ((leftSize == 2) && (rightSize == 1)) {
	  for (x=0; x<4; x++)
		 for (y=0; y<4; y++)
			for (z=0; z<4; z++) {
			   energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,z), BasePenalty(bp_pos_i+2,y), BasePenalty(bp_pos_j-1,x));
			   energy_help = Sum_MaxDouble(energy_help, interior_loop_1_2_energy[384*bp_before+96*z+24*x+4*bp_new+y]);
			   if (energy_help < min)
				  min = energy_help;
			}
	  energy = min;
   } else if ((leftSize == 2) && (rightSize == 2)) {
	  for (x1=0; x1<4; x1++)
		 for (x2=0; x2<4; x2++)
			for (y1=0; y1<4; y1++)
			   for (y2=0; y2<4; y2++) {
				  energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,x1), BasePenalty(bp_pos_i+2,y1), BasePenalty(bp_pos_j-1,x2), BasePenalty(bp_pos_j-2,y2));
				  energy_help = Sum_MaxDouble(energy_help, interior_loop_2_2_energy[1536*bp_new+256*bp_before+64*x1+16*x2+4*y1+y2]);
				  if (energy_help < min)
					 min = energy_help;
			   }
	  energy = min;
   } else {
	  // size:
	  //**********
	  if (size <= 30)
		 energy += loop_destabilizing_energies[3*(size-1)];
	  else
	  {
		 energy += loop_destabilizing_energies[3*(30-1)];
		 energy += 1.75*RT*Math.log(size/30.0);
	  }

	  // terminal mismatches a both closings:
	  //***************************************
	  // Attention: there are dependencies, if the IL has size 1 at one side of the loop,
	  //            then this base is involved in the terminal mismatches at both closings
	  if (leftSize == 1) {     //        base_i_before - base_j_before
							   //  i1                                      i2
							   //                                          i3
							   //                 bp_i - bp_j
		 for (i1=0; i1<4; i1++)
			for (i2=0; i2<4; i2++)
			   for (i3=0; i3<4; i3++) {
				  energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i1), BasePenalty(bp_pos_j-rightSize,i2), BasePenalty(bp_pos_j-1,i3));
				  energy_help = Sum_MaxDouble3(energy_help, mismatch_energies_interior[64*base_j_before+16*i2+4*base_i_before+i1], mismatch_energies_interior[64*bp_i+16*i1+4*bp_j+i3]);

				  if (energy_help < min)
					 min = energy_help;
			   }
		 energy = Sum_MaxDouble(energy, min);
	  } else if (rightSize == 1) {
							  //        base_i_before - base_j_before
							  //  i1                                      i3
							  //  i2
							  //                 bp_i - bp_j
		 for (i1=0; i1<4; i1++)
			for (i2=0; i2<4; i2++)
			   for (i3=0; i3<4; i3++) {
				  energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+leftSize,i1), BasePenalty(bp_pos_i+1,i2), BasePenalty(bp_pos_j-1,i3));
				  energy_help = Sum_MaxDouble3(energy_help, mismatch_energies_interior[64*base_j_before+16*i3+4*base_i_before+i1], mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i3]);

				  if (energy_help < min)
					 min = energy_help;
			   }
		 energy = Sum_MaxDouble(energy,min);
	  } else {
							  //        base_i_before - base_j_before
							  //  i1                                      i3
							  //  i2                                      i4
							  //                 bp_i - bp_j
		   for (i1=0; i1<4; i1++)
			  for (i3=0; i3<4; i3++) {
				 energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+leftSize,i1), BasePenalty(bp_pos_j-rightSize,i3), mismatch_energies_interior[64*base_j_before+16*i3+4*base_i_before+i1]);
				 if (energy_help < min)
					min = energy_help;
			  }
		   energy = Sum_MaxDouble(energy,min);

		   min = Number.MAX_VALUE;
		   for (i2=0; i2<4; i2++)
			  for (i4=0; i4<4; i4++) {
				 energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i2), BasePenalty(bp_pos_i-1,i4), mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i4]);
				 if (energy_help < min)
					min = energy_help;
			  }
		   energy = Sum_MaxDouble(energy,min);
	  }

	  // asymmetry-penalty:
	  //********************
	  // (for 1x2 IL not extra penalty)
	  if (leftSize != rightSize)
	  {
		 asym = 0.5 * Math.abs(leftSize-rightSize);
		 if (asym > 3.0)
			asym = 3.0;
		 energy = Sum_MaxDouble(energy,asym);
	  }
   }
   return energy;
}


/************************************************************

finds the energy of an interior loop, if the assignments of
both closing BPs and all free bases are known

(only energy fragments that are important for the difference 
of two assignments of the loop)

************************************************************/

function InteriorLoopEnergy(leftSize, rightSize, int_seq, pos_i, pos_j)
{
   var energy = 0.0;
   var x1, x2, y1, y2, i1, i2, i3, i4;

   var bp_i = int_seq[pos_i];
   var bp_j = int_seq[pos_j];
   var bp_i_before = int_seq[pos_i+leftSize+1];
   var bp_j_before = int_seq[pos_j-rightSize-1];

   var bp_new = BP2int(bp_i, bp_j);
   var bp_before = BP2int(bp_i_before, bp_j_before);

   // no bulges and loops of size 0:
   //--------------------------------------------
   if ((leftSize == 0) || (rightSize == 0)) {
	  throw("No interior loop!");
   }

   // Special cases:
   //------------------------------------------------------
   if ((leftSize == 1) && (rightSize == 1)) {
	  x1 = int_seq[pos_i+1];
	  y1 = int_seq[pos_j-1];
	  energy += interior_loop_1_1_energy[96*bp_new+24*x1+4*bp_before+y1];
   } else if ((leftSize == 1) && (rightSize == 2)) {
	  x1 = int_seq[pos_i+1];
	  y1 = int_seq[pos_j-2];
	  y2 = int_seq[pos_j-1];
	  energy += interior_loop_1_2_energy[384*bp_new+96*y1+24*x1+4*bp_before+y2];
   }
   else if ((leftSize == 2) && (rightSize == 1)) {
	  x1 = int_seq[pos_i+1];
	  x2 = int_seq[pos_i+2];
	  y1 = int_seq[pos_j-1];
	  energy += interior_loop_1_2_energy[384*bp_before+96*x1+24*y1+4*bp_new+x2];
   }
   else if ((leftSize == 2) && (rightSize == 2)) {
	  x1 = int_seq[pos_i+1];
	  x2 = int_seq[pos_i+2];
	  y1 = int_seq[pos_j-2];
	  y2 = int_seq[pos_j-1];
	  energy += interior_loop_2_2_energy[1536*bp_new+256*bp_before+64*x1+16*y2+4*x2+y1];
   } else {
	  // size: (is not important here, since this function is just used for the calculation of a difference
	  // of two assignments of the same IL)
	  //----------------------------------------
	  /*if (size <= 30)
		 energy += loop_destabilizing_energies[3*(size-1)];
	  else
	  {
		 energy += loop_destabilizing_energies[3*(30-1)];
		 energy += 1.75*RT*log((double)(size/30.0));
	  }*/

	  // terminal mismatches a both closings:
	  //----------------------------------------
	  // Attention: there are dependencies, if the IL has size 1 at one side of the loop,
	  //            then this base is involved in the terminal mismatches at both closings
	  if (leftSize == 1) {      //        bp_i_before - bp_j_before
								//  i1                                      i2
								//                                          i3
								//               bp_i - bp_j
		 i1 = int_seq[pos_i+1];
		 i2 = int_seq[pos_j-rightSize];
		 i3 = int_seq[pos_j-1];
		 energy += mismatch_energies_interior[64*bp_j_before+16*i2+4*bp_i_before+i1] + mismatch_energies_interior[64*bp_i+16*i1+4*bp_j+i3];
	  } else if (rightSize == 1) {
									 //        bp_i_before - bp_j_before
									 //  i1                                      i3
									 //  i2
									 //               bp_i - bp_j
		 i1 = int_seq[pos_i+leftSize];
		 i2 = int_seq[pos_i+1];
		 i3 = int_seq[pos_j-1];
		 energy += mismatch_energies_interior[64*bp_j_before+16*i3+4*bp_i_before+i1] + mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i3];
	  } else {
							  //        bp_i_before - bp_j_before
							  //  i1                                      i3
							  //  i2                                      i4
							  //               bp_i - bp_j
		   i1 = int_seq[pos_i+leftSize];
		   i2 = int_seq[pos_i+1];
		   i3 = int_seq[pos_j-rightSize];
		   i4 = int_seq[pos_j-1];
		   energy += mismatch_energies_interior[64*bp_j_before+16*i3+4*bp_i_before+i1] + mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i4];
	  }

	  // asymmetry-penalty: (not important here)
	  //------------------------
	  // (for 1x2 IL not extra penalty)
	  /*if (leftSize != rightSize)
	  {
		 asym = 0.5 * abs(leftSize-rightSize);
		 if (asym > 3.0)
			asym = 3.0;
		 energy+=asym;
	  }*/
   }
   return energy;
}

/*******************************************************
 identifies connections concerning the closing base pairs
 of the stems and concerning the free bases adjacent to them

 stem_ends = contains all stem ends and their positions in BP_Order, we are interested in the last stem_num
 order_pos_of_closing_bp = pos. of the closing BP in BP_Order
 BaseConnections = contains row-wise the connections concerning bases
 BasePairConnections = contains row-wise the connections concerning base pairs

 Return: vector including BPs in the ML (= predecessors of the closing BP of the ML)
*******************************************************/

function MultiLoopConnections(stem_num, stem_ends, order_pos_of_closing_bp, BaseConnections, BasePairConnections)
{
   var ml_stem_ends = new Array(stem_num);    //stem ends in the ML
   var ml_between_stems = new Array(stem_num + 1);  //contains the number of bases between the stems (in each coord. the bases "before" the
                                                                    //stem, only at the last coord. of bases_between_stems the free bases after
                                                                    //the last stem and before the closing BP are stored

   var one_connection = new Array();         //help vector for one connection

   // get just the vareresting stem ends from stem_ends and erase them there,
   // thereby their order is reversed. BPs (i,j) having a small i are first now
   for (var i=0; i<stem_num; i++)
   {
      ml_stem_ends[i] = stem_ends[stem_ends.length - 1];
      stem_ends.pop();
   }

   // number of free bases between the stem-ending base pairs
   ml_between_stems[0] = BP_Order[ml_stem_ends[0]][0] - BP_Order[order_pos_of_closing_bp][0] - 1;
   for (var i=1; i<stem_num; i++)
      ml_between_stems[i] = BP_Order[ml_stem_ends[i]][0] - BP_Order[ml_stem_ends[i-1]][1] - 1;
   ml_between_stems[stem_num] = BP_Order[order_pos_of_closing_bp][1] - BP_Order[ml_stem_ends[stem_num-1]][1] - 1;

   // finding all connections: (concerning free bases):
   //*****************************************************

   // if there's at least one free base upstream the first stem
   if (ml_between_stems[0] != 0)
      one_connection.push(BP_Order[order_pos_of_closing_bp][0]+1);

   for (var i=0; i<stem_num; i++)
   {
      // if there is only one free base between the stems and if there is at least one free base after the next stem
      if (ml_between_stems[i] == 1)
      {
         if (ml_between_stems[i+1] != 0)
            one_connection.push(BP_Order[ml_stem_ends[i]][1]+1);
      }
      else
      {
         BaseConnections.push(one_connection);
         one_connection.length = 0;
         if (ml_between_stems[i] != 0)
            one_connection.push(BP_Order[ml_stem_ends[i]][0]-1);
         if (ml_between_stems[i+1] != 0)
            one_connection.push(BP_Order[ml_stem_ends[i]][1]+1);
      }
   }

   if (ml_between_stems[stem_num] == 1)
   {
      BaseConnections.push(one_connection);
      one_connection.length = 0;
   }
   else
   {
      BaseConnections.push(one_connection);
      one_connection.length = 0;
      if (ml_between_stems[stem_num] != 0)
      {
         one_connection.push(BP_Order[order_pos_of_closing_bp][1]-1);
      }
      // one_connection might be empty, but this is done to ensure that BaseConnections and BasePairConnections
      // are always of the same size
      BaseConnections.push(one_connection);
      one_connection.length = 0;
   }

   // finding all connections: (concerning stems):
   //***********************************************

   one_connection.push(order_pos_of_closing_bp);

   for (var i=0; i<stem_num; i++)
   {
      if (ml_between_stems[i] == 1)
         one_connection.push(ml_stem_ends[i]);
      else
      {
         BasePairConnections.push(one_connection);
         one_connection.length = 0;
         one_connection.push(ml_stem_ends[i]);
      }
   }

   if (ml_between_stems[stem_num] == 1)
   {
      one_connection.push(order_pos_of_closing_bp);
      BasePairConnections.push(one_connection);
      one_connection.length = 0;
   }
   else
   {
      BasePairConnections.push(one_connection);
      one_connection.length = 0;
      one_connection.push(order_pos_of_closing_bp);
      BasePairConnections.push(one_connection);
      one_connection.length = 0;
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

function MultiLoopConnections2(stem_ends, order_pos_of_closing_bp, BaseConnections, BaseSizes, BasePairConnections, BasePairSizes)
{
   var i;
   var stem_num = BP_Order[order_pos_of_closing_bp][3];
   var num_bp_con = 0;   //current number of BasePairConnections
   var num_base_con = 0; //current number of BaseConnections
   var bp, base;         //current position in the current connection

   var ml_stem_ends = new Array(stem_num);  //stem ends of the ML
   var ml_between_stems = new Array(stem_num + 1); //contains the number of bases between the stems (in each coord. the bases "before" the
                                                                   //stem, only at the last coord. of bases_between_stems the free bases after
                                                                   //the last stem and before the closing BP are stored

   // get just the vareresting stem ends from stem_ends and erase them there,
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

function BestConnectionEnergy(base_connection, pair_connection, bp_at_pair_connection)
{
   var energy = 0.0;
   var min = Number.MAX_VALUE;
   var base_energy;

   for (var i=0; i< base_connection.length; i++)
   {
      //for each base, test, whether it is adjacent to one or two stems, then: minimize energy

      //bases in base_connection[i] may be neighbored to base pairs in
      //pair_connection[i-1] [i] or/and [i+1] => testing (and also whether the respecting field in the vector exists)
      min = Number.MAX_VALUE;

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }
      }

      if (i < pair_connection.length)
      {
         // if the left stem is the closingBP
         if (BP_Order[pair_connection[i]][0]+1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located left of (= before) the stem
         if (BP_Order[pair_connection[i]][0]-1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i]][1]+1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i][1]+4*bp_at_pair_connection[i][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // if the right stem is the closingBP
         if (BP_Order[pair_connection[i]][1]-1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

      }

      if ( i+1 < pair_connection.length)
      {
         // usual case: base is located left of (= before) the stem
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i+1][1]+4*bp_at_pair_connection[i+1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
                  min = base_energy;
            }

         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i+1]][1]-1 == base_connection[i])
            for (var j=0; j<4; j++)
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

function ConnectionEnergy(base_connection, base_size, pair_connection, bp_size, int_seq)
{
	var energy = 0.0;
	var i;
	var min;

	for(i = 0; i < base_size; i++)
	{
    	//for each base, test, whether it is adjacent to one or two stems, then: minimize energy
	
    	// bases in base_connections[i] may be neighbored to base pairs in
    	// pair_connection[i-1] [i] or/and [i+1] => testing (and also whether the respecting field in the vector exists)
    	min = Number.MAX_VALUE;

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

function MLBestEnergy(bp_pos, stem_ends)
{
   var min = Number.MAX_VALUE;
   var MLenergy, energy;
   var offset = 3.4;
   var free_base_penalty = 0.0;
   var helix_penalty = 0.4;
   var energy_help = 0.0, energy_unbound, energy_stems;
   var bp_i, bp_j;
   var pair_size, base_size;
   var i;
   var ML_vorgaenger; // BP_Order pos. of the stems in the ML
   var BaseConnections;     // all connections conc. free bases
   var BasePairConnections; // all connections conc. stems

   var bp_at_pair_connection;      // base pair assignment for a connection
   var MIN_bp_at_pair_connection;  // base pair assignment for a connection with min. free energy

   stem_ends.push(bp_pos-1);    //previous BP is a stem ending BP in the ML
   //BP_Order[bp_pos][3]+1, since also closingBP-stems counts
   MLenergy = offset + free_base_penalty * BP_Order[bp_pos][2] + helix_penalty * (BP_Order[bp_pos][3]+1);
   ML_vorgaenger = MultiLoopConnections(BP_Order[bp_pos][3], stem_ends, bp_pos, BaseConnections, BasePairConnections);

   // since the closing BP of the ML can be included in two connections, its assignment is fixed before testing the assignments of the other BPs
   for (var bp_assign=0; bp_assign<6; bp_assign++)
   {
      //store predecessors in Trace, assignment still unknown
      for (var vg=0; vg<BP_Order[bp_pos][3]; vg++)
         Trace[bp_pos][bp_assign][vg][0] = ML_vorgaenger[vg];

      energy = MLenergy;  // energy

      for (i=0; i<BaseConnections.length; i++) //BaseConnections.size() == BasePairConnections.size()
      {
         pair_size = BasePairConnections[i].length;
         base_size = BaseConnections[i].length;

         bp_at_pair_connection = new Array(pair_size);
         MIN_bp_at_pair_connection = new Array(pair_size);
         for (var a=0; a<pair_size; a++)
         {
            bp_at_pair_connection[a] = new Array(2);
            MIN_bp_at_pair_connection[a] = new Array(2);
            for (var b=0; b<2; b++)
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

         var start_pos; // if the closing BP is the first one in the connection, start_pos for assigning bp_at_pair_connection
                        // is not 0, but 1
         var start_power;  // if the closing BP is the last in the connection, start_power for assigning bp_at_pair_connection
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

         min = Number.MAX_VALUE;

         //this is the case, if the connection only consists of the closing BP, then the following for-loop is not done.
         //Nevertheless, an adjacent free base has to be analyzed.
         if ((pair_size == 1) && (start_power == pair_size-3))
         {
            if (base_size == 0)
               min = 0;
            else
            {
               // energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here)
               energy_unbound = BestConnectionEnergy(BaseConnections[i], BasePairConnections[i], bp_at_pair_connection);
               min = energy_unbound;
            }
         }

         for (var int_anz = 0; int_anz < Math.pow(6, start_power+1); int_anz++)
         {
            var power = start_power;
            var value = int_anz;
            var pos = start_pos;
            var bp_help;

            // create an assignment
            while (power >= 0)
            {
               bp_help = value/ (Math.pow(6,power));

               BP2_2(bp_help, bp_i, bp_j);
               bp_at_pair_connection[pos][0]=bp_i;
               bp_at_pair_connection[pos][1]=bp_j;

               value = value % Math.pow(6, power);

               power--;
               pos++;
            }

            // energy of a connection with fixed bp-assignments of the stems (only the energy of the free bases is added here, energy of 
            // the stem will be added later
            energy_unbound = BestConnectionEnergy(BaseConnections[i], BasePairConnections[i], bp_at_pair_connection);

            // add the energy of the stems (depending on the the closing BP), BUT: not closing BP!
            energy_stems = 0.0;
            for (var j=0; j<pair_size; j++)
               if (BasePairConnections[i][j] != bp_pos)
                  energy_stems = Sum_MaxDouble(energy_stems,D[BasePairConnections[i][j]][BP2int(bp_at_pair_connection[j][0],bp_at_pair_connection[j][1])]);

            energy_help = Sum_MaxDouble(energy_stems, energy_unbound);

            if (min > energy_help)
            {
               min = energy_help;
               for (var a=0; a<pair_size; a++)
               {
                   MIN_bp_at_pair_connection[a][0] = bp_at_pair_connection[a][0];
                   MIN_bp_at_pair_connection[a][1] = bp_at_pair_connection[a][1];
               }
            }
         }// for int_anz

         // if no valid assignment was found, set a random one (this is the case, if in the loop a BP exists, that has no valid assignment
         // concerning the constraints --> usually this error is filtered before)
         if (min == Number.MAX_VALUE)
            for (var a=0; a<pair_size; a++)
            {
               var rand_bp = RandomBasePair();
               BP2_2(rand_bp, MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);
            }

         // add minimal energy that arises because of the free bases in the ML
         energy = Sum_MaxDouble(energy,min);

         // store best assignment of the predecessors in the traceback
         for (var a=0; a<pair_size; a++)
            if (BasePairConnections[i][a] != bp_pos)
               for (var vg = 0; vg<BP_Order[bp_pos][3]; vg++)
                  if (BasePairConnections[i][a] == Trace[bp_pos][bp_assign][vg][0])
                     Trace[bp_pos][bp_assign][vg][1] = BP2int(MIN_bp_at_pair_connection[a][0], MIN_bp_at_pair_connection[a][1]);

      } // for i
      
      BP2_2(bp_assign, bp_i, bp_j);
      D[bp_pos][bp_assign] = Sum_MaxDouble3(energy, Zero_or_StemEndAU(bp_pos, bp_assign),PairPenalty(bp_pos, bp_i, bp_j));

   } //for bp_assign

   for (i=0; i<BaseConnections.length; i++)
      BaseConnections[i].length = 0;
   BaseConnections.length = 0;

   for (i=0; i<BasePairConnections.length; i++)
      BasePairConnections[i].length = 0;
   BasePairConnections.length = 0;
}

/****************************************************************************************

finds energy of the ML (for a given sequence!!) (only the fraction that is relevant
for the difference of two assignments of it)

****************************************************************************************/

function MLEnergy(bp_pos, int_seq)
{
   var energy = 0.0;
   //double offset = 3.4;
   //double free_base_penalty = 0.0;
   //double helix_penalty = 0.4;
   var energy_unbound;
   var pair_size, base_size;
   var i, j;
   var BaseConnections;     // all connections, concerning free bases
   var BasePairConnections; // all connections, concerning stems
   var BaseSizes;            // sizes of the base connections
   var BasePairSizes;        // sizes of the stem connections

   var bp_assign; // assignment of the closing BP
   var num_conns;
   var stem_ends;
   var con_size = BP_Order[bp_pos][3]+2; // in order that base and basepairconn. have the same size, no conn. can be bigger than #stems+2 and all together
                                                  // there can not be more (since only adjacent free base are taken into account)
   stem_ends = FindStemEnds(bp_pos);

   BaseConnections = new Array(con_size);
   BaseSizes = new Array(con_size);
   for (i=0; i<con_size; i++)
   {
      BaseConnections[i] = new Array(con_size);
      BaseSizes[i] = -1;
      for (j=0;j<con_size; j++)
         BaseConnections[i][j] = -1;
   }

   BasePairConnections = new Array(con_size);
   BasePairSizes = new Array(con_size);
   for (i=0; i<con_size; i++)
   {
      BasePairConnections[i] = new Array(con_size);
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

   for (i=0; i< num_conns; i++) //BaseConnections.size() == BasePairConnections.size()
   {
      pair_size = BasePairSizes[i];
      base_size = BaseSizes[i];

      // energy of one connection (just energy of the free bases)
      energy_unbound = ConnectionEnergy(BaseConnections[i], base_size, BasePairConnections[i], pair_size, int_seq);
      energy = Sum_MaxDouble(energy, energy_unbound);
   } // for i

   //base pair penalty only for the closingBP
   energy = Sum_MaxDouble3(energy, PairPenalty(bp_pos,int_seq[BP_Order[bp_pos][0]],int_seq[BP_Order[bp_pos][1]]), Zero_or_StemEndAU(bp_pos, bp_assign));

   return energy;
}

/*****************************************************************

Finds all BPs that are stems ends in a ML

*****************************************************************/

function FindStemEnds(closingBP)
{
   var i, bp_j;
   var stem_ends;
   stem_ends = new Array(BP_Order[closingBP][3]);

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

function ConnectionBestFreeBases(base_connection, pair_connection, bp_at_pair_connection)
{
   var min = Number.MAX_VALUE;
   var base_energy;
   var base_size = base_connection.length;
   var min_base=0;
   var min_bases;
   min_bases = new Array(base_size);
   for (var i=0; i<base_size; i++)
      min_bases[i] = -1;

   for (var i=0; i<base_size; i++)
   {
      //for each base, test, whether it is adjacent to one or two stems, then: minimize energy

      //bases in base_connection[i] may be neighbored to base pairs in
      //pair_connection[i-1] [i] or/and [i+1] => testing (and also whether the respecting field in the vector exists)
      min = Number.MAX_VALUE;

      if (i-1 >= 0)  /*case BP_Order[pair_connection[i-1]][0] not possible*/
      {
         // usual case: base is located right of (= after) the stem
         if (BP_Order[pair_connection[i-1]][1]+1 == base_connection[i])
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[16*bp_at_pair_connection[i-1][1]+4*bp_at_pair_connection[i-1][0]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }
      }

      if (i < pair_connection.length)
      {
         // if the left stem is the closing BP
         if (BP_Order[pair_connection[i]][0]+1 == base_connection[i])
            for (var j=0; j<4; j++)
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
            for (var j=0; j<4; j++)
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
            for (var j=0; j<4; j++)
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
            for (var j=0; j<4; j++)
            {
               base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_at_pair_connection[i][0]+4*bp_at_pair_connection[i][1]+j],BasePenalty(base_connection[i],j));
               if (base_energy < min)
               {
                  min = base_energy;
                  min_base = j;
               }
            }
      }

      if ( i+1 < pair_connection.length)
      {
         // usual case: base is located left of (= before) the stem
         if (BP_Order[pair_connection[i+1]][0]-1 == base_connection[i])
            for (var j=0; j<4; j++)
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
            for (var j=0; j<4; j++)
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