#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;


int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
	OptStruct oiflag = {"i", 0, NULL, NULL, "flag for inverse binarization"}; options.push_back(&oiflag);	
	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	ParStruct pvalue = {"threshold", NULL, "threshold"}; parameters.push_back(&pvalue);
	
	
	if (!parsecmdline("src_binarize", "binarize image values", argc, argv, options, parameters))
		return 0;
	
	
	
	
	
	//! input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	
	//! process
	libUSTG::cflimage output = input.binarize(atof(pvalue.value), oiflag.flag);

	
	//! save
	output.save(pout.value);
	

}

