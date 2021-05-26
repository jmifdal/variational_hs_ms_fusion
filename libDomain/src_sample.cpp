#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char **argv)
{


	
	vector <OptStruct *> options; 
	OptStruct oa = {"a", 0, NULL, NULL,  "flag for aglomeration"};  options.push_back(&oa);	
	OptStruct og = {"g:", 0, NULL, NULL, "flag for gaussian convolution and std"};  options.push_back(&og);	
	
    
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	ParStruct psample = {"s", NULL, "integer sampling step"}; parameters.push_back(&psample);
	
	
	if (!parsecmdline("src_sample", "Samples image of factor s", argc, argv, options, parameters))
		return 0;
		
	
	int sample_factor =  atoi(psample.value);
		
	
	
	// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	libUSTG::cflimage output;
	if (og.flag)
	{
		float std = atof(og.value);
        input.subSampleConv(sample_factor,  BOUNDARY_CONDITION_SYMMETRIC, std, &output);
	
	} else  if (oa.flag)
	{	
		input.subSampleAglomeration(sample_factor, &output);
		
	} else 
	{
		input.subSample(sample_factor, &output);
	}

	
	output.save( pout.value); 
	return 1;
	

}
