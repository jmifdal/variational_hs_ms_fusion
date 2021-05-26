#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;


int main(int argc, char **argv)
{

	
	vector <OptStruct *> options;
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);

	
	if (!parsecmdline("src_convert", "convert image format", argc, argv, options, parameters))
		return 0;


	
	//! Read input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	input.save(pout.value); 
	return 1;
	
}
