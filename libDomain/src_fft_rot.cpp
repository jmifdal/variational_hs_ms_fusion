#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char* argv[])
{
	
	
	
	vector <OptStruct *> options;
	OptStruct of  = {"f", 0, NULL,   NULL, "if activated same origin and width, height used"}; options.push_back(&of);
	OptStruct os  = {"s", 0, NULL,   NULL, "if activated image simmetrized before translation"}; options.push_back(&os);
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout   = {"out",   NULL, "output image"}; parameters.push_back(&pout);
	ParStruct pa   = {"a",   NULL, "angle in degrees"}; parameters.push_back(&pa);
	ParStruct px   = {"x",   NULL, "translation in x"}; parameters.push_back(&px);
	ParStruct py   = {"y",   NULL, "translation in y"}; parameters.push_back(&py);
	
	
	if (!parsecmdline("src_fft_rot", "Rotates image", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);

	//! Parameters	
	float angle =  atof(pa.value);
	float tx =  atof(px.value);
	float ty =  atof(py.value);
	
	
	
	//! Process
	input.fftRot(angle, tx, ty, of.flag, os.flag);

	
	//! Save
	input.save( pout.value); 
	
}


