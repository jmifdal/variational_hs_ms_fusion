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
	ParStruct pinput = {"image1", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image2", NULL, "image2"}; parameters.push_back(&pinput2);
	ParStruct pmask = {"mask", NULL, "a value"}; parameters.push_back(&pmask);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	if (!parsecmdline("src_linear_combination_mask", "Computes mask*I1 + (1-mask)*I2",
                      argc, argv, options, parameters))
		return 0;
	
	
	
	
	// Input
	libUSTG::cflimage input;
	input.load(pinput.value);

	
	// Input
	libUSTG::cflimage input2;
	input2.load(pinput2.value);
	
	
	if (!input.isSameSize(input2)) { printf("images of different size or number of channels"); }
	
    
    libUSTG::cflimage mask;
    mask.load(pmask.value);
    
	
    // Process
	libUSTG::cflimage output = input;
	for (int ii = 0; ii < input.whc(); ii++)
	{
		output[ii] = mask[ii % input.wh()] * input[ii] + (1.0f - mask[ii % input.wh()])* input2[ii];
	}
	
 	output.save(pout.value); 
	return 1;
	
	
}
