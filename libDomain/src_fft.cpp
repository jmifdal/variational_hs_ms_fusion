#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;

int main(int argc, char* argv[])
{



	vector <OptStruct *> options;
	OptStruct oR = {"r:", 0, NULL, NULL, "real part"};  options.push_back(&oR);
    OptStruct oI= {"i:", 0, NULL, NULL, "imaginary part"};  options.push_back(&oI);
	OptStruct oM = {"m:", 0, NULL, NULL, "modul"};  options.push_back(&oM);
    OptStruct oL = {"l:", 0, NULL, NULL, "log of modulus"};  options.push_back(&oL);
    
    
    
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	
	
	if (!parsecmdline("src_fft", "Computes Fourier transform", argc, argv, options, parameters))
		return 0;

	

	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

    if (input.c() > 1)
    {
        printf("warning :: image with more than one channels, using only first one\n");
    }
    
    
    libUSTG::flimage iim(input.w(), input.h()); iim = 0.0f;
    libUSTG::flimage ore(input.w(), input.h()); ore = 0.0f;
    libUSTG::flimage oim(input.w(), input.h()); oim = 0.0f;
    libUSTG::flimage modul(input.w(), input.h()); modul = 0.0f;
    
    
	
    //! Compute fft
    int i_flag = 0;
    libUSTG::fft2d(input.v(), iim.v(), ore.v(), oim.v(), i_flag, input.w(), input.h());
    
    
    //! Compute modul
    for (int ii=0; ii < input.wh(); ii++)  modul[ii] = sqrtf( ore[ii]*ore[ii] + oim[ii]*oim[ii]);

    
    
    

	
	//! Save
    if (oR.flag) ore.save(oR.value);
    if (oI.flag) oim.save(oI.value);
    
    if (oM.flag) modul.save(oM.value);
    if (oL.flag)
    {
        for (int ii=0; ii < input.wh(); ii++)  modul[ii] = (float) log(1. + (double) modul[ii]);
        modul.save(oL.value);
    }
    
    
}


