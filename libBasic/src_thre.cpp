#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;
int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
	OptStruct om = {"m:", 0, NULL, NULL, "minimum value"}; options.push_back(&om);	
	OptStruct oM = {"M:", 0, NULL, NULL, "maximum value"}; options.push_back(&oM);	

	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_thre", "Shrinks image values", argc, argv, options, parameters))
		return 0;
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	
	//! Parameters
	float fminim = 0.0f;
	if (om.flag) fminim = atof(om.value);
	else fminim = input.min();

	
	float fmaxim = 0.0f;
	if (oM.flag) fmaxim = atof(oM.value);
	else fmaxim = input.max();

	
	
	//! Process
	input.thre(fminim, fmaxim);	
	
	input.save(pout.value);
	
}

