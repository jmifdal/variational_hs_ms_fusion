#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;

/// Difference between images
int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image1", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image2", NULL, "image"}; parameters.push_back(&pinput2);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_add", "Add the values of two images per pixel", argc, argv, options, parameters))
		return 0;



	
    

	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	libUSTG::cflimage input2;
	input2.load(pinput2.value);

	 

	if (! input.isSameSize(input2))
	{
		printf("images of different size or number of channels\n");
		return 0;
	}

	
	
	//! Process
	libUSTG::cflimage output = input;
	for(int i=0; i < input.whc(); i++)  output[i] += input2[i];

    
	output.save( pout.value); 
	return 1;
	
}
