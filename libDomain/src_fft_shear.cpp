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
	OptStruct oD  = {"d:", 0, "0",   NULL, "orientation (0: horizontal, 1:vertical"}; options.push_back(&oD);
	OptStruct oS  = {"s", 0, NULL,   NULL, "if activated symmetrization applied"}; options.push_back(&oS);
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout   = {"out",   NULL, "output image"}; parameters.push_back(&pout);
	ParStruct pa  = {"angle",   NULL, "shear factor"}; parameters.push_back(&pa);
	ParStruct px  = {"translation",  NULL, "translation factor"}; parameters.push_back(&px);
	
	
	if (!parsecmdline("src_fft_shear", "Shears image", argc, argv, options, parameters))
		return 0;
	

	
	//! Parameters	
	int flagNoCenter =  of.flag;
	int flagSymmetric = oS.flag;
	int orientation = atoi(oD.value);
	if (orientation != 0) orientation= 1;
	
	
	float shear =  atof(pa.value);
	float translation = atof(px.value);
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	//! Process		
	input.fftShear(shear, orientation, translation, flagNoCenter, flagSymmetric);
	
	
	//! Save
	input.save( pout.value); 

}


