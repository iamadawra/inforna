/***********************************************************
Dependencies Begin Here!
***********************************************************/

/**********************************************************
 translates a BP (given as an integer) to two integerBases
**********************************************************/
function BP2_2(bp, base_i, base_j)
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
function BP2_2(bp, base_i, base_j)
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

function Sum_MaxDouble3(sum1, sum2, sum3)
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
	Vvar bp_at_pair_connection;

	var base_energy, energy_help;

	var min_dang;
	var min_free_base;
	if (min_result[0] == MAX_DOUBLE){
		bp_assign = RandomBasePair();
	}else{
		bp_assign = min_result[0];
	}
	BP2_2(bp_assign,bp_assign_i,bp_assign_j);

	var vorg;
	vorg = new Array(BP_Order[numBP][3]);
	for (var vg=0; vg< Maximum(BP_Order[numBP][3],1); vg++){
		vorg[vg] = Trace[numBP][bp_assign][vg][0];
		BP2_2(Trace[numBP][bp_assign][vg][1],bp_i,bp_j);
		int_seq[BP_Order[vorg[vg]][0]] = bp_i;
		int_seq[BP_Order[vorg[vg]][1]] = bp_j;
	}

	if (BP_Order[numBP-1][0] > 0){
		min_dang = MAX_DOUBLE;
		min_free_base = -1;
		for (var free_base = 0; free_base < 4; free_base++){
			base_energy = Sum_MaxDouble(single_base_stacking_energy[64+16*bp_assign_j+4*bp_assign_i+free_base], BasePenalty(BP_Order[numBP-1][0]-1,free_base));
			if (min_dang > base_energy){
				min_dang = base_energy;
				min_free_base = free_base;
			}
		}

		int_seq[BP_Order[numBP-1][0]-1] = min_free_base;
		if (BP_Order[numBP-1][0] > 1){
			for (var vor = 0; vor < BP_Order[numBP][0]-1; vor++){
				int_seq[vor] = SetFreeBase(vor);
			}
		}
	}

	if (BP_Order[numBP-1][1] < brackets.length - 1){
		EndConnections(BP_Order[numBP][3], vorg, BaseConnections, BasePairConnections);

		for (i=0; i < BaseConnections.size(); i++){
			pair_size = BasePairConnections[i].size();

			bp_at_pair_connection = new Array(pair_size);
			for (var a=0; a<pair_size; a++){
				bp_at_pair_connection[a] = new Array(2);
				for(var b=0; b<2; b++){
					bp_at_pair_connection[a][b] = 0;
				}
			}

			for (var bp=0; bp<pair_size; bp++){
				bp_at_pair_connection[bp][0] = int_seq[BP_Order[BasePairConnections[i][bp]][0]];
				bp_at_pair_connection[bp][1] = int_seq[BP_Order[BasePairConnections[i][bp]][1]];
			}

			var best_end_connection_bases;
			best_end_connection_bases = EndConnectionBestFreeBases(BaseConnections[i], BasePairConnections[i], bp_at_pair_connection);
			
			base_size = BaseConnections[i].size();

			for(var pos=0; pos<base_size; pos++){
				int_seq[BaseConnections[i][pos]] = best_end_connection_bases[pos];
			}

		}

		for (i=0; i < BaseConnections.size(); i++){
			BaseConnections[i].clear();
		}
		BaseConnections.clear();

		for (i=0; i<BasePairConnections.size(); i++){
			BasePairConnections[i].clear();
		}
		BasePairConnections.clear();
	}

	for (var bp_pos=numBP-1; bp_pos>-1; bp_pos--){
	/*Traceback is done stepwise, base pair per base pair. Usually, the previous BP in BP_Order
	is the predecessor of the current one, except for MLs! MLs have more than one predecessors.
	This doesn't matter, the BPs are processed in the order of BP_Order, nevertheless. 
	
	all assignments of base pairs that are already fixed by the traceback, are written in int_seq,
	i.e. concerning MLs: the closingBP and the last BPs of the stems. If the traceback (running
	through BP_Order) reaches these BPs, they are already fixed in int_seq. For go on in the traceback,
	choose the fixed assignment (given in int_seq) and find the right assignment of the
	predecessor(s) with the help of Trace.*/

	/*************************************************************************/

	bp_assign = BP2int(int_seq[BP_Order[bp_pos][0]],int_seq[BP_Order[bp_pos][1]]);
	hairpin_loop = false;

	var bp_pos_i = BP_Order[bp_pos][0];
	var bp_pos_j = BP_Order[bp_pos][1];

	//finding the assignments of the predecessors and store in int_seq
	for (var vg=0; vg<Maximum(BP_Order[bp_pos][3],1); vg++)
	{
 		if (Trace[bp_pos][bp_assign][vg][1] != -1) //i.e. it has a predecessor (no closing BP of a HL)
 		{
    		BP2_2(Trace[bp_pos][bp_assign][vg][1], bp_i, bp_j);
    		int_seq[BP_Order[Trace[bp_pos][bp_assign][vg][0]][0]] = bp_i;
    		int_seq[BP_Order[Trace[bp_pos][bp_assign][vg][0]][1]] = bp_j;
 		}
 		else //a HL follows
 		{
    		hairpin_loop = true;
    		break;
 		}
	}

    valid_hairpin_loop = false;
    BP2_2(bp_assign,bp_assign_i,bp_assign_j);

    //*****************************************************************************
    //*****************************************************************************
    //assign HAIRPINLOOP
    //*************************
    if (hairpin_loop)
    {
 		//find loop size:
 		hairpin_loop_size = BP_Order[bp_pos][1] - BP_Order[bp_pos][0] - 1;
 		if (hairpin_loop_size == 3)
 		{
    		// in case of a triloop, the energy only depends on the size (=3), exceptions: CCC => +1.4; GGG => -2.2
    		// thus firstly, triloop are assigned randomly, maybe considering all possibilities later
    		while (valid_hairpin_loop == false)
    		{
       			//finding an assignment
       			for (var pos=1; pos<=3; pos++)
          			int_seq[bp_pos_i+pos] = SetFreeBase(bp_pos_i+pos);

       			//test, whether the assignment gives the energy that is stored in D (done for making the CCC-loop valid)
       			//penalties for AU-closingBP and GU-closingBP and CCC-Loop, negative penalty for GGG-Loop (GGG missed out, since also done in Vienna Package)
       			energy = loop_destabilizing_energies[3*hairpin_loop_size-1];

       			// not energy fragment for the closing BP, but penalties for it and the bases in the loop
       			energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,int_seq[bp_pos_i+1]), BasePenalty(bp_pos_i+2,int_seq[bp_pos_i+2]), BasePenalty(bp_pos_i+3,int_seq[bp_pos_i+3]), PairPenalty(bp_pos, bp_assign_i, bp_assign_j));

       			energy = Sum_MaxDouble(energy_help,energy);

       			if ((bp_assign == 0) || (bp_assign == 3) || (bp_assign == 4) || (bp_assign == 5))
          			energy = Sum_MaxDouble(energy,terminalAU);

       			if ((int_seq[bp_pos_i+1] == 1) && (int_seq[bp_pos_i+2] == 1) && (int_seq[bp_pos_i+3] == 1))
          			energy = Sum_MaxDouble(energy,Ctriloop);
       			//else if ((int_seq[BP_Order[bp_pos][0]+1] == 2) && (int_seq[BP_Order[bp_pos][0]+2] == 2) && (int_seq[BP_Order[bp_pos][0]+3] == 2))
       			//   energy = Sum_MaxDouble(energy,Gtriloop);

       			//terminal mismatch energy, if the HLclosing BP is the last one in a stem
       			if ((StackEnd(bp_pos)) && ((bp_assign == 0) || (bp_assign == 3) || (bp_assign == 4) || (bp_assign == 5)))
          			energy = Sum_MaxDouble(energy,terminalAU);

       			//test, whether HL is valid
       			if (D[bp_pos][bp_assign] == energy)
          			valid_hairpin_loop = true;
       			// dealing with GGG triloops, missed out, since also done in Vienna Package
		        /*else if (Sum_MaxDouble(D[bp_pos][bp_assign], Gtriloop) == energy)
		        {
		           int_seq[BP_Order[bp_pos][0]+1] = 2;
		           int_seq[BP_Order[bp_pos][0]+2] = 2;
		           int_seq[BP_Order[bp_pos][0]+3] = 2;
		           valid_hairpin_loop = true;
		        }*/
    		}
 		}
 		else if (hairpin_loop_size == 4)
 		{
		    //allocate memory and assign BP
		    var tetra_plus_closing;
		    tetra_plus_closing = new Array(6);
		    tetra_plus_closing[0] = int2char(bp_assign_i);
		    tetra_plus_closing[5] = int2char(bp_assign_j);

		    energy = loop_destabilizing_energies[3*hairpin_loop_size-1];

		    // consider all cases of the tetraloop (4*4*4*4) and add term-mismatch and
		    // possibly a bonus, furthermore add the penalties for all 4 bases
		    // (even if there could be some equally good loop assignments, choose one randomly)
		    //      i2   j2
		    //   i           j
		    //     bp_i-bp_j
		    //       .   .

		    min = MAX_DOUBLE;
		    for (var i=0; i<4; i++)
		       for (var j=0; j<4; j++)
		          for (var i2=0; i2<4; i2++)
		             for (var j2=0; j2<4; j2++)
		             {
		                energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,i), BasePenalty(bp_pos_i+2,i2), BasePenalty(bp_pos_j-2,j2),BasePenalty(bp_pos_j-1,j));

		                tetra_plus_closing[1] = int2char(i);
		                tetra_plus_closing[2] = int2char(i2);
		                tetra_plus_closing[3] = int2char(j2);
		                tetra_plus_closing[4] = int2char(j);

		                energy_help = Sum_MaxDouble3(energy_help, tetra_loop_energy(tetra_plus_closing), mismatch_energies_hairpin[64*bp_assign_i+16*i+4*bp_assign_j+j]);

		                // If the HPLoop would only consist of C's, a penalty of 0.3*size+1.6 has to be added.
		                // This means that just in case that all bases in the loop are restricted to C by the
		                // constraints the energy has to be corrected.
		                // NOT TAKEN INTO ACCOUNT, SINCE THIS IS NOT TAKEN INTO ACCOUNT IN THE VIENNA PACKAGE AS WELL.
		                //if ((i==1) && (i2==1) && (j2==1) && (j==1))
		                //   energy_help = Sum_MaxDouble(energy_help,0.3*hairpin_loop_size+1.6);

		                if (energy_help < min)
		                {
		                   min = energy_help;
		                   min_i = i;
		                   min_i2 = i2;
		                   min_j2 = j2;
		                   min_j = j;
		                }
		             }

		    energy = Sum_MaxDouble(energy,min);
		    int_seq[bp_pos_i+1] = min_i;
		    int_seq[bp_pos_i+2] = min_i2;
		    int_seq[bp_pos_j-2] = min_j2;
		    int_seq[bp_pos_j-1] = min_j;
		    free(tetra_plus_closing);

		    //terminal mismatch energy, if the HLclosing BP is the last one in a stem
		    if ((StackEnd(bp_pos)) && ((bp_assign == 0) || (bp_assign == 3) || (bp_assign == 4) || (bp_assign == 5)))
		       energy = Sum_MaxDouble(energy,terminalAU);
		}
 		else //>4
		{
    		while (valid_hairpin_loop == false)
    		{
       			energy = loop_destabilizing_energies[3*Minimum(hairpin_loop_size,30)-1];
       			if (hairpin_loop_size > 30)
          			energy += 1.75*RT*log(hairpin_loop_size/30.0);

       			min = MAX_DOUBLE;
       			for (var i=0; i<4; i++)
          			for (var j=0; j<4; j++)
          			{
             			// no energy fragment by the closing BP, but penalties have to be taken into account
             			energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,i), BasePenalty(bp_pos_j-1,j), PairPenalty(bp_pos, bp_assign_i, bp_assign_j), mismatch_energies_hairpin[64*bp_assign_i+16*i+4*bp_assign_j+j]);

             			if (energy_help < min)
             			{
                			min = energy_help;
                			min_i = i;
                			min_j = j;
             			}
          			}

       			energy = Sum_MaxDouble(energy,min);
       			int_seq[bp_pos_i+1] = min_i;
       			int_seq[bp_pos_j-1] = min_j;

       			//find an assignment for the loop (the first and the last base are fixed by min. terminal mismatch,
       			//thus, starting at the second pos. end at the last but one pos.
       			for (var pos=2; pos<=hairpin_loop_size-1; pos++)
          			int_seq[bp_pos_i+pos] = SetFreeBase(bp_pos_i+pos);

       			// If the HPLoop would only consist of C's, a penalty of 0.3*size+1.6 has to be added.
       			// NOT TAKEN INTO ACCOUNT, SINCE THIS IS NOT TAKEN INTO ACCOUNT IN THE VIENNA PACKAGE AS WELL.

       			// bool onlyCs = true;
       			// for (int p=bp_pos_i+1; p<bp_pos_j; p++)
       			//    if ( int_seq[p] != 1) 
       			//    {
       			//       onlyCs = false;
       			//       break;
       			//    }

       			// if (onlyCs)
       			//    energy = Sum_MaxDouble(energy,0.3*hairpin_loop_size+1.6);

       			//terminal mismatch energy, if the HLclosing BP is the last one in a stem
       			if ((StackEnd(bp_pos)) && ((bp_assign == 0) || (bp_assign == 3) || (bp_assign == 4) || (bp_assign == 5)))
          			energy = Sum_MaxDouble(energy,terminalAU);

       			if (D[bp_pos][bp_assign] == energy)
          			valid_hairpin_loop = true;
    		}
 		}
	}

      //**********************************************************************************
      //**********************************************************************************
      //MULTILOOP
      //*****************
      //fixed the assignments of the free bases in the ML
      else if (BP_Order[bp_pos][3] > 1)
      {
         var ML_vorgaenger; // pos. in BP_Order of the BPs in the ML, = predecessors of the closing BP

         // finding all stem ends of the ML (stored in the information in Trace)
         vector <int> stem_ends;
         for (var vg=BP_Order[bp_pos][3]-1; vg>=0; vg--)
            stem_ends.push_back(Trace[bp_pos][bp_assign][vg][0]);

         ML_vorgaenger = MultiLoopConnections(BP_Order[bp_pos][3], stem_ends, bp_pos, BaseConnections, BasePairConnections);

         for (i=0; i<BaseConnections.size(); i++) //BaseConnections.size() == BasePairConnections.size()
         {
            pair_size = BasePairConnections[i].size();

            bp_at_pair_connection = (int**) malloc(sizeof(int*)*pair_size);
            for (var a=0; a<pair_size; a++)
            {
               bp_at_pair_connection[a] = (int*) malloc(sizeof(int)*2);
               for (var b=0; b<2; b++)
                  bp_at_pair_connection[a][b] = 0;
            }

            /* for bp_at_pair_connection choose the assignment given in Trace*/
            for (var bp=0; bp<pair_size; bp++)
            {
               if (BasePairConnections[i][bp] == bp_pos)
               {
                  bp_at_pair_connection[bp][0] = bp_assign_i;
                  bp_at_pair_connection[bp][1] = bp_assign_j;
               }
               else
               {
                  bp_at_pair_connection[bp][0] = int_seq[BP_Order[BasePairConnections[i][bp]][0]];
                  bp_at_pair_connection[bp][1] = int_seq[BP_Order[BasePairConnections[i][bp]][1]];
               }
            }

            var best_connection_bases;
            best_connection_bases = ConnectionBestFreeBases(BaseConnections[i], BasePairConnections[i], (const int**) bp_at_pair_connection);

            base_size = BaseConnections[i].size();
            for (var pos=0; pos<base_size; pos++)
               int_seq[BaseConnections[i][pos]] = best_connection_bases[pos];

            for (var c=0; c<pair_size; c++)
               free(bp_at_pair_connection[c]);

         } // for i

         for (i=0; i<BaseConnections.size(); i++)
            BaseConnections[i].clear();
         BaseConnections.clear();

         for (i=0; i<BasePairConnections.size(); i++)
            BasePairConnections[i].clear();
         BasePairConnections.clear();
      }

      //*******************************************************************************************
      //*******************************************************************************************
      //remaining structural elements including free bases (BULGE, INTERIOR LOOP)
      //************************************************************
      else if (bp_pos > 0) //then there are predecessors and BL or IL can arise
      {
         left_loop_size = BP_Order[bp_pos-1][0] - BP_Order[bp_pos][0] - 1;
         right_loop_size = BP_Order[bp_pos][1] - BP_Order[bp_pos-1][1] - 1;

         //STACK (do nothing, since there are no free bases)
         //*************
         if ((left_loop_size == 0) && (right_loop_size == 0))
         {
         }

         //LEFT BULGE
         //*************
         else if ((left_loop_size != 0) && (right_loop_size == 0))
         {
            //whatever assignments of the bases in the BL ==> Constraints checked in SetFreeBase
            for (var pos=1; pos<=left_loop_size; pos++)
               int_seq[BP_Order[bp_pos][0]+pos] = SetFreeBase(BP_Order[bp_pos][0]+pos);
         }

         // RIGHT BULGE
         //****************
         else if ((left_loop_size == 0) && (right_loop_size != 0))
         {
            //whatever assignments of the bases in the BL ==> Constraints checked in SetFreeBase
            for (int pos=1; pos<=right_loop_size; pos++)
               int_seq[BP_Order[bp_pos-1][1]+pos] = SetFreeBase(BP_Order[bp_pos-1][1]+pos);
         }


         // INTERIOR LOOP
         //****************
         else
         {
            var bp_before = Trace[bp_pos][bp_assign][0][1];
            var bp_before_i, bp_before_j;
            BP2_2(bp_before,bp_before_i,bp_before_j);
            min = MAX_DOUBLE;

            // special cases:
            //****************************************************
            if ((left_loop_size == 1) && (right_loop_size == 1))
            {
               for (var x=0; x<4; x++)
                  for (var y=0; y<4; y++)
                  {
                     energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,x), BasePenalty(bp_pos_j-1,y), interior_loop_1_1_energy[96*bp_assign+24*x+4*bp_before+y]);
                     if (energy_help < min)
                     {
                        min = energy_help;
                        min1 = x;
                        min2 = y;
                     }
                  }
               int_seq[bp_pos_i+1] = min1;
               int_seq[bp_pos_j-1] = min2;
            }
            else if ((left_loop_size == 1) && (right_loop_size == 2)) /*x(i-1)-x(i)=2 and y(i)-y(i-1)=3*/
            {
               for (var x=0; x<4; x++)
                  for (var y=0; y<4; y++)
                     for (var z=0; z<4; z++)
                     {
                        energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,x), BasePenalty(bp_pos_j-1,y), BasePenalty(bp_pos_j-2,z), interior_loop_1_2_energy[384*bp_assign+96*z+24*x+4*bp_before+y]);
                        if (energy_help < min)
                        {
                           min = energy_help;
                           min1 = x;
                           min2 = y;
                           min3 = z;
                        }
                     }
               int_seq[bp_pos_i+1] = min1;
               int_seq[bp_pos_j-1] = min2;
               int_seq[bp_pos_j-2] = min3;
            }
            else if ((left_loop_size == 2) && (right_loop_size == 1)) /*x(i-1)-x(i)=3 and y(i)-y(i-1)=2*/
            {
               for (var x=0; x<4; x++)
                  for (var y=0; y<4; y++)
                     for (var z=0; z<4; z++)
                     {
                        energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,z), BasePenalty(bp_pos_i+2,y), BasePenalty(bp_pos_j-1,x), interior_loop_1_2_energy[384*bp_before+96*z+24*x+4*bp_assign+y]);
                        if (energy_help < min)
                        {
                           min = energy_help;
                           min1 = x;
                           min2 = y;
                           min3 = z;
                        }
                     }
               int_seq[bp_pos_i+1] = min3;
               int_seq[bp_pos_i+2] = min2;
               int_seq[bp_pos_j-1] = min1;
            }
            else if ((left_loop_size == 2) && (right_loop_size == 2))
            {
               for (var x1=0; x1<4; x1++)
                  for (var x2=0; x2<4; x2++)
                     for (var y1=0; y1<4; y1++)
                        for (var y2=0; y2<4; y2++)
                        {
                           energy_help = Sum_MaxDouble4(BasePenalty(bp_pos_i+1,x1), BasePenalty(bp_pos_i+2,y1), BasePenalty(bp_pos_j-1,x2), BasePenalty(bp_pos_j-2,y2));
                           energy_help = Sum_MaxDouble(energy_help, interior_loop_2_2_energy[1536*bp_assign+256*bp_before+64*x1+16*x2+4*y1+y2]);
                           if (energy_help < min)
                           {
                              min = energy_help;
                              min1 = x1;
                              min2 = x2;
                              min3 = y1;
                              min4 = y2;
                           }
                        }
               int_seq[bp_pos_i+1] = min1;
               int_seq[bp_pos_i+2] = min3;
               int_seq[bp_pos_j-1] = min2;
               int_seq[bp_pos_j-2] = min4;
            }

            //from now on: just interactions with closingBPs, no special cases
            //**********************************************************************
            else
            {
               // terminal mismatches at both closings:
               //***************************************
               // Attention: there are dependencies, if the IL has size 1 at one side of the loop,
               //            then this base is involved in the terminal mismatches at both closings
               if (left_loop_size == 1)      //          bp_before_i - bp_before_j
               {                       //  i1                                      i2
                                       //                                          i3
                                       //                bp_assign_i - bp_assign_j
                  for (var i1=0; i1<4; i1++)
                     for (var i2=0; i2<4; i2++)
                       for (var i3=0; i3<4; i3++)
                       {
                          energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i1), BasePenalty(bp_pos_j-right_loop_size,i2), BasePenalty(bp_pos_j-1,i3));

                          energy_help = Sum_MaxDouble3(energy_help, mismatch_energies_interior[64*bp_before_j+16*i2+4*bp_before_i+i1],
                          mismatch_energies_interior[64*bp_assign_i+16*i1+4*bp_assign_j+i3]);

                          if (energy_help < min)
                          {
                             min = energy_help;
                             min1 = i1;
                             min2 = i2;
                             min3 = i3;
                          }
                       }

                  int_seq[bp_pos_i+1] = min1;
                  int_seq[bp_pos_j-right_loop_size] = min2;
                  int_seq[bp_pos_j-1] = min3;

                   //remaining bases are assigned randomly
                  for (var i=2; i<=right_loop_size-1; i++)
                     int_seq[bp_pos_j-i] = SetFreeBase(bp_pos_j-i);
               }

               else if (right_loop_size == 1)//          bp_before_i - bp_before_j
               {                       //  i1                                      i3
                                       //  i2
                                       //                bp_assign_i - bp_assign_j
                  for (var i1=0; i1<4; i1++)
                     for (var i2=0; i2<4; i2++)
                        for (var i3=0; i3<4; i3++)
                        {
                           energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+left_loop_size,i1), BasePenalty(bp_pos_i+1,i2), BasePenalty(bp_pos_j-1,i3));

                           energy_help = Sum_MaxDouble3(energy_help, mismatch_energies_interior[64*bp_before_j+16*i3+4*bp_before_i+i1], mismatch_energies_interior[64*bp_assign_i+16*i2+4*bp_assign_j+i3]);

                           if (energy_help < min)
                           {
                              min = energy_help;
                              min1 = i1;
                              min2 = i2;
                              min3 = i3;
                           }
                        }

                  int_seq[bp_pos_j-1] = min3;
                  int_seq[bp_pos_i+1] = min2;
                  int_seq[bp_pos_i+left_loop_size] = min1;

                   //remaining bases are assigned randomly
                  for (var i=2; i<=left_loop_size-1; i++)
                     int_seq[bp_pos_i+i] = SetFreeBase(bp_pos_i+i);
               }

               else                    //        bp_before_i - bp_before_j
               {                       //  i1                                      i3
                                       //  i2                                      i4
                                       //        bp_assign_i - bp_assign_j
                  for (var i1=0; i1<4; i1++)
                     for ( var i3=0; i3<4; i3++)
                     {
                        energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+left_loop_size,i1), BasePenalty(bp_pos_j-right_loop_size,i3), mismatch_energies_interior[64*bp_before_j+16*i3+4*bp_before_i+i1]);
                        if (energy_help < min)
                        {
                           min = energy_help;
                           min1 = i1;
                           min3 = i3;
                        }
                     }
                  int_seq[BP_Order[bp_pos-1][0]-1] = min1;
                  int_seq[BP_Order[bp_pos-1][1]+1] = min3;

                  min = MAX_DOUBLE;
                  for (var i2=0; i2<4; i2++)
                     for ( var i4=0; i4<4; i4++)
                     {
                        energy_help = Sum_MaxDouble3(BasePenalty(bp_pos_i+1,i2), BasePenalty(bp_pos_j-1,i4), mismatch_energies_interior[64*bp_assign_i+16*i2+4*bp_assign_j+i4]);
                        if (energy_help < min)
                        {
                           min = energy_help;
                           min2 = i2;
                           min4 = i4;
                        }
                     }
                  int_seq[bp_pos_i+1] = min2;
                  int_seq[bp_pos_j-1] = min4;

                  //remaining bases are assigned randomly
                  for (var i=2; i<=left_loop_size-1; i++)
                     int_seq[bp_pos_i+i] = SetFreeBase(bp_pos_i+i);
                  for (var i=2; i<=right_loop_size-1; i++)
                     int_seq[bp_pos_j-i] = SetFreeBase(bp_pos_j-i);
               }
            }
         }
      }
   } //for bp_pos
   return int_seq;
}
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


