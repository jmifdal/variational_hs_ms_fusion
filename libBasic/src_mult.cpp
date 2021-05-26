#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;


int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
	OptStruct oa = {"a:", 0, NULL , NULL,"a value"};  options.push_back(&oa);	
	OptStruct oA = {"A:", 0, NULL , NULL,"image value"};  options.push_back(&oA);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput1 = {"image1", NULL, "image"}; parameters.push_back(&pinput1);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);


	if (!parsecmdline("src_mult", "Multiplies by a value or another image", argc, argv, options, parameters))
		return 0;

	
	if (!oa.flag && !oA.flag)
	{
		printf("error :: at least one  of options -a or -A must be selected"); exit(-1);
	}

	
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput1.value);
	
	int width = input.w();
	int height = input.h();
	int nchannels = input.c();
	
	
	libUSTG::cflimage input2;
	if (oA.flag) {

		input2.load(oA.value);
		int width2 = input2.w();
		int height2 = input2.h();
		int nchannels2 = input2.c();

		if (width != width2 || height != height2 || nchannels != nchannels2) 
		{
		
			printf("error :: images of different size"); exit(-1);
		}
		
	}
	
	

	libUSTG::cflimage output(width,height,nchannels);
	if (oA.flag)
	{
	
		for(int i = 0; i < nchannels; i++)
		{
			float *u = input.v(i);
			float *v = input2.v(i);
			float *o = output.v(i);
			
			for (int j=0; j < width*height; j++)
			{
				o[j] = u[j] * v[j];
			}
		}
			
	} else {
		
		float a = atof(oa.value);
		for(int i = 0; i < nchannels; i++)
		{
			float *u = input.v(i);
			float *o = output.v(i);
			
			for (int j=0; j < width*height; j++)
			{
				o[j] = u[j] * a;
			}
		}
		
	}


	

	output.save( pout.value); 
	return 1;
}
