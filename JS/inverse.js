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








function getOrder(structure){
	var bp_pos; 
	var next_pos;
	var k=0;
	var stack_len;
	var help_num = new Array(5);
	var results_ML = new Array(2);
	var bpTable;
	bpTable = getPairmap(structure);

	// ElementStruct = Element_Structure(right_Element_Structure(structure));
	// <--- Convert this line of CPP code. No library functionality

	BP_Order = new Array(numBP + 1); // Still no love for numBP

	for (var i=0; i<numBP+1; i++){
		BP_Order[i] = new Array(4);
		for (var j=0; j<4; j++){
			BP_Order[i][j] = 0;
		}
	}

	//Finding the BP Order
	//***************************
	bp_pos = numBP - 1;
	for (var i=0; i<structure.length; i++){
		if ((bpTable[i]!= -1) && (bpTable[i]>i)){
			if (bp_pos < 0){
				throw("Problem to get the order!");
				return false; //In place of Exiting the function
			}
			BP_Order[bp_pos][0] = i;
			BP_Order[bp_pos][1] = bpTable[i];
			bp_pos--;
		}
	}

	//Finding closing BPs of MLs
	//****************************
	bp_pos = numBP - 1;
	var i = 0;
	while(i<ElementStruct.length){
		if ((ElementStruct[i]=="(") && (ElementStruct[i+1]=="S")){
			k=0;
			stack = true;
			while (ElementStruct[i+k+2] != '('){ 
				help_num[k] = ElementStruct[i+k+2];
				k++;
			}
			help_num[k] = '\0';
			stack_len = parseInt(help_num);
			bp_pos = stack_len;
		}

		next_pos = i+k+2+1;

		if ((ElementStruct[next_pos]=='M') && (stack == true)){
			results_ML = StemsInML(next_pos-1, bp_pos);

			bp_pos = results_ML[1];
			i = results_ML[0];
		}
		stack = false;
		i++;
	}
	var stems_and_freeBases;

	//stems_and_freeBases = ClosingStructure(bpTable);<--- Convert this line of CPP code. No library functionality

	BP_Order[numBP][2] = stems_and_freeBases[0];
	BP_Order[numBP][3] = stems_and_freeBases[1];
}

function StemsInML(first_pos, BP_start_pos){
	var k = 0;
	var pos = first_pos;
	var bp_pos = BP_start_pos;
	var multi = true;
	var help_num = new Array(5);
	var stack_len = 0;
	var stems = 0;
	var pos_return;
	pos_return = new Array(2);

	pos += 2;

	k = 0;
	while (ElementStruct[pos+k] != '('){
		help_num[k] = ElementStruct[pos+k];
		k++;
	}
	help_num[k] ='\0';
	BP_Order[BP_start_pos+1][2] = parseInt(help_num);

	while (ElementStruct[pos] != '('){
		pos++;
	}

	while (multi){
		if ((ElementStruct[pos] == '(') && (ElementStruct[pos+1]=='S')){
			k = 0;
			while (ElementStruct[pos+k+2] != '('){
				help_num[k] = ElementStruct[pos+k+2];
				k++;
			}
			help_num[k] = '\0';
			stack_len += parseInt(help_num);
			bp_pos -= parseInt(help_num);
			pos += k+2;
		}

		else if ((EleemntStruct[pos]=='S') && (ElementStruct[pos-1]!= '(')){
			k = 0;
			while(ElementStruct[pos+k+1] != ')'){
				help_num[k] = ElementStruct[pos+k+1];
				k++;
			}
			help_num[k] ='\0';
			stack_len -= parseInt(help_num);
			if (stack_len == 0){
				stems++;
				pos += k+2;
			}

		}

		else if ((ElementStruct[pos]=='(') && (ElementStruct[pos+1]=='M')){
			pos_return = StemsInML(pos,bp_pos);
			pos = pos_return[0] + 1;
			bp_pos = pos_return[1];
		} 

		else if ((ElementStruct[pos]=='M') && (ElementStruct[pos-1]!='(')){
			multi = false;
			while (ElementStruct[pos] != ')'){
				pos++;
			}
		}

		else{
			pos++;
		}

		BP_Order[BP_start_pos+1][3] = stems;
		pos_return[0] = pos;
		pos_return[1] = bp_pos;

		return pos_return;
	}
}

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



function Recursion(){
	var anz_vorgaenger;
	var min_vorgaenger;
	var best_int_seq;
	var energy_help;

	D = new Array(numBP+1);
	Trace = new Array(numBP+1);

	for(var i=0; i<numBP;i++){
		D[i] = new Array(6);
		Trace[i] = new Array(6);

		anz_vorgaenger = max(1, BP_Order[i][3]);

		for (var j = 0; j<6; j++;){
			D[i][j] = 0.0;
			Trace[i][j] = new Array(anz_vorgaenger);
			for (var k=0; k<anz_vorgaenger; k++){
				Trace[i][j][k] = new Array(2);
				for (var u=0; u<2; u++){
					Trace[i][j][k][u] = -1;
				}
			}	
		}
	}

	var new_stem = true;
	var min = Number.MAX_VALUE;
	var loop_size, left_loop_size, right_loop_size;
	var bp_i, bp_j;
	//VECTOR LINE - Line 255
	var bp_pos;
	for(bp_pos=0; bp_pos < numBP; bp_pos++){
		if (new_stem){
			if (bp_pos){
				loop_size = BP_Order[bp_pos][1] - BP_Order[bp_pos][0]-1;
				for (var bp_assign=0; bp_assign < 6; bp_assign++){
					BP2_2(bp_assign, bp_i, bp_j);
					D[bp_pos][bp_assign] = Sum_MaxDouble3()
				}
			}
		}
	}
}

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
