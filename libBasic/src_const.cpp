#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;


int main(int argc, char **argv)
{
	
	
	
	vector <OptStruct *> options;
	OptStruct ow = {"w:", 0,  "256", NULL, "image width"};  options.push_back(&ow);	
	OptStruct oh = {"h:", 0,  "256", NULL, "image height"};  options.push_back(&oh);	
	OptStruct on = {"n:", 0,  "1", NULL, "number of channels"};  options.push_back(&on);
	OptStruct ov = {"v:", 0,  "0", NULL, "plane value"};  options.push_back(&ov);
	
	
	vector<ParStruct *> parameters;
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_const", "Creates constant image", argc, argv, options, parameters))
		return 0;
	
	
	//! Parameters
	int height = atoi(oh.value);
	int width = atoi(ow.value);
	int nchannels = atoi(on.value);
	float fValue = atof(ov.value);

	
	//! Process
	libUSTG::cflimage output(width, height, nchannels);

	output = fValue;
	output.save(pout.value);
	
}

