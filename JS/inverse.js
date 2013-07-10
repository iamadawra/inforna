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
      throw("Wrong BasepairNr.!");
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
					D[bp_pos][bp_assign] = Sum_MaxDouble3(BestHairpinLoopEnergy(bp_pos,loop_size, bp_i, bp_j), Zero_or_StemEndAU(bp_pos, bp_assign), PairPenalty(bp_pos, bp_i, bp_j));

				}
			}
			else{
				if (BP_Order[bp_pos][1] < BP_Order[bp_pos-1][0]){
					//VECTOR PUSH BACK FUNCTIONALITY
					loop_size = BP_Order[bp_pos][1] - BP_Order[bp_pos][0]-1;
					for (var bp_assign=0; bp_assign<6;bp_assign++){
						BP2_2(bp_assign, bp_i, bp_j);
						D[bp_pos][bp_assign] = Sum_MaxDouble3(BestHairpinLoopEnergy(bp_pos, loop_size, bp_i, bp_j), Zero_or_StemEndAU(bp_pos, bp_assign), PairPenalty(bp_pos, bp_i, bp_j));
					}
				}
				else{
					MLBestEnergy(bp_pos, stem_ends);
				}
			}
			new_stem = false;
		}
		else{
			left_loop_size = BP_Order[bp_pos-1][0] - BP_Order[bp_pos][0] - 1;
			right_loop_size = BP_Order[bp_pos][1] - BP_Order[bp_pos-1][1] - 1;

			for (var bp_assign=0; bp_assign<6; bp_assign++){
				BP2_2(bp_assign, bp_i, bp_j);
				min = MAX_DOUBLE;
				min_vorgaenger = -1;

				// Stack
				if ((left_loop_size == 0) && (right_loop_size == 0)){
					for (var bp_before=0; bp_before < 6; bp_before++){
						energy_help = Sum_MaxDouble(StackingEnergy(bp_i, bp_j, bp_before),D[bp_pos-1][bp_before]);
						if (min > energy_help){
						 min = energy_help;
						 min_vorgaenger = bp_before;
						}
				   }
					D[bp_pos][bp_assign] = Sum_MaxDouble3(min, Zero_or_StemEndAU(bp_pos, bp_assign), PairPenalty(bp_pos, bp_i, bp_j));
					Trace[bp_pos][bp_assign][0][0] = bp_pos-1;
					//if all assignments of the predecessor are set to MAX_DOUBLE, choose an assignment of the predecessor randomly
					if (min == MAX_DOUBLE)
				    	Trace[bp_pos][bp_assign][0][1] = RandomBasePair();
					else
				    	Trace[bp_pos][bp_assign][0][1] = min_vorgaenger;
				}

				// left bulge
				else if ((left_loop_size != 0) && (right_loop_size == 0)){
					for (var bp_before=0; bp_before < 6; bp_before++){
				    	energy_help = Sum_MaxDouble(BulgeEnergy(left_loop_size, bp_i, bp_j, bp_before), D[bp_pos-1][bp_before]);
				    	if (min > energy_help){
				    		min = energy_help;
				    		min_vorgaenger = bp_before;
				    	}
					}
					D[bp_pos][bp_assign] = Sum_MaxDouble3(min, Zero_or_StemEndAU(bp_pos, bp_assign), PairPenalty(bp_pos, bp_i, bp_j));
					Trace[bp_pos][bp_assign][0][0] = bp_pos-1;
					//if all assignments of the predecessor are set to MAX_DOUBLE, choose an assignment of the predecessor randomly
					if (min == MAX_DOUBLE)
				    	Trace[bp_pos][bp_assign][0][1] = RandomBasePair();
					else
				    	Trace[bp_pos][bp_assign][0][1] = min_vorgaenger;
				}

				// right bulge
				else if ((left_loop_size == 0) && (right_loop_size != 0)){
					for (var bp_before=0; bp_before < 6; bp_before++){
				    	energy_help = Sum_MaxDouble(BulgeEnergy(right_loop_size, bp_i, bp_j, bp_before), D[bp_pos-1][bp_before]);
				    	if (min > energy_help){
				        	min = energy_help;
				        	min_vorgaenger = bp_before;
				    	}
					}
					D[bp_pos][bp_assign] = Sum_MaxDouble3(min, Zero_or_StemEndAU(bp_pos, bp_assign), PairPenalty(bp_pos, bp_i, bp_j));
					Trace[bp_pos][bp_assign][0][0] = bp_pos-1;
					//if all assignments of the predecessor are set to MAX_DOUBLE, choose an assignment of the predecessor randomly
					if (min == MAX_DOUBLE)
				    	Trace[bp_pos][bp_assign][0][1] = RandomBasePair();
					else
				    	Trace[bp_pos][bp_assign][0][1] = min_vorgaenger;
				}

				// interior loop
				else{
				   for (var bp_before=0; bp_before < 6; bp_before++){
				    	energy_help = Sum_MaxDouble(BestInteriorLoopEnergy(bp_pos, left_loop_size, right_loop_size, bp_i, bp_j, bp_before),D[bp_pos-1][bp_before]);
				    	if (min > energy_help){
				        	min = energy_help;
				        	min_vorgaenger = bp_before;
				    	}
					}
					D[bp_pos][bp_assign] = Sum_MaxDouble3(min, Zero_or_StemEndAU(bp_pos, bp_assign), PairPenalty(bp_pos, bp_i, bp_j));
					Trace[bp_pos][bp_assign][0][0] = bp_pos-1;
					//if all assignments of the predecessor are set to MAX_DOUBLE, choose an assignment of the predecessor randomly
					if (min == MAX_DOUBLE)
				    	Trace[bp_pos][bp_assign][0][1] = RandomBasePair();
					else
				    	Trace[bp_pos][bp_assign][0][1] = min_vorgaenger;
				}
			} // end for bp_assign
		} // end if else stem_end
		if (StackEnd(bp_pos)){
			new_stem = true;
		}
	}

	externBestEnergy();

	best_int_seq = Traceback();

	for (var i=0; i<struct_len; i++){
		if (best_int_seq[i] == -1){
			best_int_seq[i] = SetFreeBase(i);
		}
	}

	best_char_seq = new Array(struct_len+1);
	for (var i=0; i<struct_len; i++){
		best_char_seq[i] = int2char(best_int_seq[i]);
	}
	best_char_seq[struct_len] = '\0'; // IS THIS LINE NECESSARY??? Array controller basically...
	var min_result = MiniVec(D[numBP], 6); //MiniVec ????????????????????

	return min_result[1];
}


/****************************************************************************
*****************************************************************************

 Doing the traceback:

 *****************************************************************************
 ****************************************************************************/


function Traceback(){
	var seq_len = brackets.length;
	var int_seq;
	int_seq = new Array(seq_len);
	for (var i=0; i<seq_len; i++){
		int_seq[i] = -1;
	}
	var int_result;
	min_result = MiniVec(D[numBP], 6); //MINIVEC ?????????
	var bp_assign, bp_assign_i, bp_assign_j, min_i=0, min_j=0, min_i2=0, min_j2=0, min1=0, min2=0, min3=0, min4=0;
	var bp_i, bp_j;
	var hairpin loop;
	var valid_hairpin_loop;
	var left_loop_size, right_loop_size, hairpin_loop_size;
	var energy, min=0;

	var pair_size, base_size;
	var i;
	//Vector declarations 


	//Milestone July 9 - 8.15pm
}







/*********************************************************/
/* Random initialization                                 */
/*********************************************************/

function Random_Init(){
	var min_en;
	var bp_i, bp_j;
	var sum_constraints, rand_base, rand_pair;

	best_char_seq = new Array(struct_len+1);
	best_char_seq[struct_len] = '\0';  // AGAIN. This \0 is annoying ..

	var bpTable;  // stores for each pos. the bound pos. in the BP (or -1 if unbound)
	bpTable = make_BasePair_Table(brackets);

	for (var i=0; i<struct_len; i++){
		if (bpTable[i] == -1){
	     
	    	sum_constraints = SumVec(seq_constraints[i], 4);
	    	
	    	rand_base = RandomBase(sum_constraints) + 1;
	     
	    	var ones = 0;
	    	var column = -1;
	    	while (ones < rand_base){
	    		column++;
	        	if (seq_constraints[i][column] == 1)
	        		ones++;
	     	}
	    	//column is the randomly chosen base assignment (0=A, 1=C, 2=G, 3=U)
	    	best_char_seq[i] = int2char(column);         
		}
	    else if (bpTable[i] > i){
	        var bp_constraints; //stores for all bp-assignments whether they are allowed or not
	        bp_constraints = (int*) malloc(sizeof(int)*6);
	        for (int bp=0; bp<6; bp++)
	        	bp_constraints[bp] = 0;
	        // finding allowed pairs
	        if ((seq_constraints[i][0] == 1) && (seq_constraints[bpTable[i]][3] == 1))
	        	bp_constraints[0] = 1;
	        if ((seq_constraints[i][1] == 1) && (seq_constraints[bpTable[i]][2] == 1))
	        	bp_constraints[1] = 1;
	        if ((seq_constraints[i][2] == 1) && (seq_constraints[bpTable[i]][1] == 1))
	        	bp_constraints[2] = 1;
	        if ((seq_constraints[i][3] == 1) && (seq_constraints[bpTable[i]][0] == 1))
	        	bp_constraints[3] = 1;
	        if ((seq_constraints[i][2] == 1) && (seq_constraints[bpTable[i]][3] == 1))
	        	bp_constraints[4] = 1;
	        if ((seq_constraints[i][3] == 1) && (seq_constraints[bpTable[i]][2] == 1))
	        	bp_constraints[5] = 1;
	            
	        sum_constraints = SumVec(bp_constraints, 6);
	         
	        rand_pair = RandomBasePair(sum_constraints) + 1;
	         
	        var ones = 0;
	        var column = -1;
	        while (ones < rand_pair){
	        	column++;
	            if (bp_constraints[column] == 1)
	            	ones++;
	        }
	        //column is the randomly chosen base pair assignment (0=AU, 1=CG, 2=GC, 3=UA, 4=GU, 5=UG)
	        BP2_2(column, bp_i, bp_j);
	        best_char_seq[i] = int2char(bp_i);
	        best_char_seq[bpTable[i]] = int2char(bp_j);
	    }
	}
	   
	var test_str, string;
	test_str = new Array(struct_len+1);
	string = new Array(struct_len+1);
	strcpy(string, best_char_seq); //STRCOPY IN JS ???????
	min_en = fold(string, test_str);
	return min_en;
}


