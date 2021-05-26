#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;


int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
	OptStruct om = {"m:", 0, NULL, NULL, "r value"};  options.push_back(&om);	
	OptStruct oM = {"M:", 0, NULL, NULL, "r value"};  options.push_back(&oM);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image1", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	if (!parsecmdline("src_float_2_char", "converts float to char image", argc, argv, options, parameters))
		return 0;
	
	
	
	
	// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	// Parameters
	float minim;
	if (om.flag)   
		minim = atof(om.value);
	else
		minim = input.min();
	
	
	float maxim;
	if (oM.flag)   
		maxim = atof(oM.value);
	else
		maxim = input.max();
	
	
	// Process
	float vAmp =  maxim - minim;
	for(int i=0; i < input.whc(); i++)
	{
		
		float dif = (input[i] - minim) / vAmp;
		input[i] = 255.0 * MAX(0.0, MIN(1.0, dif));
		
		
	}
	
	
 	input.save(pout.value); 
	return 1;
	
	
}
