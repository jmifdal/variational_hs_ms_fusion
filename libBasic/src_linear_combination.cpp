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
	ParStruct pa = {"a", NULL, "a value"}; parameters.push_back(&pa);
	ParStruct pb = {"b", NULL, "b value"}; parameters.push_back(&pb);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	if (!parsecmdline("src_linear_combination", "computes a*I1 + b*I2", argc, argv, options, parameters))
		return 0;
	
	
	
	
	// Input
	libUSTG::cflimage input;
	input.load(pinput.value);

	
	// Input
	libUSTG::cflimage input2;
	input2.load(pinput2.value);
	
	
	if (!input.isSameSize(input2)) { printf("images of different size or number of channels"); }
	
	
	// Parameters
	float fA = atof(pa.value);
	float fB = atof(pb.value);
	
	
	// Process
	libUSTG::cflimage output = input;
	for (int ii=0; ii < input.whc(); ii++)
	{
		output[ii] = fA * input[ii] + fB * input2[ii];
	}
	
 	output.save(pout.value); 
	return 1;
	
	
}
