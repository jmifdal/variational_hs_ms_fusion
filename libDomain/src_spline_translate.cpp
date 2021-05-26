#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;

int main(int argc, char* argv[])
{



	vector <OptStruct *> options;
	OptStruct oB = {"b:", 0, "0",   NULL, "boundary value if not interpolated"}; options.push_back(&oB);
	OptStruct oV = {"v", 0, NULL,   NULL, "interpolate boundary value"}; options.push_back(&oV);
	

	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
    ParStruct px   = {"x",   NULL, "translation in x"}; parameters.push_back(&px);
    ParStruct py   = {"y",   NULL, "translation in y"}; parameters.push_back(&py);
    
	
	if (!parsecmdline("src_spline_translate", "Translates image by splines interpolation",
                      argc, argv, options, parameters))
		return 0;

	

	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	


	//! Process
	libUSTG::cflimage output = input;


    libUSTG::flimage transx(input.w(), input.h()); transx = -atof(px.value);
    libUSTG::flimage transy(input.w(), input.h()); transy = -atof(py.value);
    
    for (int ii = 0; ii < input.c(); ii++)
    {
        libUSTG::bicubic_interpolation_warp(input.v(ii), transx.v(), transy.v(), input.w(), input.h(), oV.flag, atof(oB.value), output.v(ii));
    }

	
	//! Save
	output.save( pout.value); 

}


