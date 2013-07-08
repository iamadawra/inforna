function getOrder(structure){
	var bp_pos; 
	var next_pos;
	var k=0;
	var stack_len;
	var help_num;
	var results_ML;
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