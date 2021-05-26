#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include"../library/libImage.h"
#include"../library/libBasic.h"



using namespace std;

int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
    OptStruct ob = {"b:", 0, "0", NULL,  "0: neumann   1: symmetry"};  options.push_back(&ob);

	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	ParStruct psigma = {"sigma", NULL, "kernel standard deviation"}; parameters.push_back(&psigma);

	
	if (!parsecmdline("src_filter_convol_gaussian", "Convolves with a gaussian kernel", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//!  Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
    //! Parameters
 	float sigma = atof(psigma.value);
	
	
	//! Process & output
	libUSTG::cflimage output;
    input.convolveGauss(sigma, atoi(ob.value), &output);
	

	//! Save output
	output.save( pout.value); 

}
