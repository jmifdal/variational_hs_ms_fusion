#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;

int main(int argc, char* argv[])
{



	vector <OptStruct *> options;
	OptStruct oZ = {"z:", 0, "1.0", NULL, "zoom x factor"};  options.push_back(&oZ);	
	OptStruct oB = {"b:", 0, "0",   NULL, "boundary value if not interpolated"}; options.push_back(&oB);
	OptStruct oV = {"v", 0, "0",   NULL, "interpolate boundary value"}; options.push_back(&oV);
	OptStruct oO = {"o:", 0, "3",   NULL, "interpolation order"}; options.push_back(&oO);
	

	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_spline_zoom", "Zoomes by Splines", argc, argv, options, parameters))
		return 0;

	

	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
    
    int order = atoi(oO.value);

    //! Process
	libUSTG::cflimage output;
    if (order != 3)
        input.upSampleGenericSplines(atof(oZ.value), order, atof(oB.value), &output);
    else
        input.upSampleCubicSplines(atof(oZ.value), oV.flag, atof(oB.value), &output);

	
	//! Save
	output.save( pout.value); 

}


