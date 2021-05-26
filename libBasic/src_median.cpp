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
	
	
	if (!parsecmdline("src_median", "computes median image value", argc, argv, options, parameters))
		return 0;
	


	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	//float fSum = 0.0f;
	for (int ii = 0; ii < input.c(); ii++)
	{
	
		float fMean = libUSTG::fpMedian(input.v(ii), input.wh());
		
		//printf("%d: mean = %f \n", ii, fMean );
        printf("%f \n",fMean );
	}

	
	
	return 1;
}
