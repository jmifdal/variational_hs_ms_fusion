#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"output", NULL, "output mirror image"}; parameters.push_back(&pout);
	ParStruct pO = {"orientation", NULL, "mirror orientation (0: horizontal  1: vertical)"}; parameters.push_back(&pO);

	if(!parsecmdline("src_mirror", "Mirror Image", argc, argv, options, parameters))
		return EXIT_FAILURE;
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	//! Process
	libUSTG::cflimage output = input.mirror(atoi(pO.value));
	
	
	//! Saving output
	output.save( pout.value); 

    
    return EXIT_SUCCESS;
}