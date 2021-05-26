#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;

int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
    OptStruct ob = {"b:", 0, "0", NULL,  "0: neumann   1: symmetry"};  options.push_back(&ob);
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	ParStruct pkernel = {"kernel", NULL, "kernel image"}; parameters.push_back(&pkernel);

	
	if (!parsecmdline("src_filter_convol", "Convolves with a kernel", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	//! Parameters
	libUSTG::flimage kernel;
	kernel.load(pkernel.value);
	
	
	//! Output
	libUSTG::cflimage output;
    input.convolve(kernel, atoi(ob.value), &output);
		
    
    
	output.save( pout.value);

}