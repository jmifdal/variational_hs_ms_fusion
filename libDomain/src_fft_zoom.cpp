#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;

int main(int argc, char* argv[])
{



	vector <OptStruct *> options;
	OptStruct ox = {"x:", 0, "1.0", NULL, "zoom x factor"};  options.push_back(&ox);	
	OptStruct oy = {"y:", 0, "1.0", NULL, "zoom y factor"};  options.push_back(&oy);	
	//OptStruct os = {"s", 0, NULL,   NULL, "if activated image simmetrized before zoom"}; options.push_back(&os);
	

	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_fft_zoom", "Zoomes by Fourier Transform", argc, argv, options, parameters))
		return 0;

	

	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	


	//! Process
	libUSTG::cflimage output;
		
	output = input.UpSampleFFT(atof(ox.value), atof(oy.value));

	
    
	//! Save
	output.save( pout.value); 

}


