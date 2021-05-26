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
	OptStruct oaflag = {"a", 0, NULL, NULL, "flag for absolute value computation"};  options.push_back(&oaflag);	
	OptStruct oM = {"m:", 0, NULL, NULL, "mask image"};  options.push_back(&oM);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image1", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image2", NULL, "image"}; parameters.push_back(&pinput2);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_diff", "computes image difference", argc, argv, options, parameters))
		return 0;



	int aflag = oaflag.flag;
	

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

	
	//! Mask
	libUSTG::flimage imask;
	float *fpMask = NULL;
	if (oM.flag) 
	{
		imask.load(oM.value); fpMask = imask.v();
	}

	
	
	//! Process
	libUSTG::cflimage output(input.w(),input.h(),input.c());
	output=0.0f;
	 
	for(int i=0; i < input.c(); i++)
	{
		
		float *igray =  input.v(i);
		float *igray2 =  input2.v(i);
		float *ogray =  output.v(i);
		
		for(int j=0; j < input.wh(); j++) 
			if (!fpMask || fpMask[j] > 0.0f) ogray[j] = igray[j] - igray2[j];
		
		if (aflag)
			for(int j=0; j < input.wh(); j++) ogray[j] = fabsf(ogray[j]);

	}
	
	
	output.save( pout.value); 
	return 1;
	
}
