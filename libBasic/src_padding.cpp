#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"output", NULL, "output image"}; parameters.push_back(&pout);
	ParStruct pw = {"w", NULL, "width of extracted zone"}; parameters.push_back(&pw);
	ParStruct ph = {"h", NULL, "height of extracted zone"}; parameters.push_back(&ph);
	ParStruct pa = {"a", NULL, "padding value"}; parameters.push_back(&pa);
	
	if(!parsecmdline("src_padding", "Pads image with value a", argc, argv, options, parameters))
		return EXIT_FAILURE;
	

	//! Parameters
	int w = atoi(pw.value);
	int h = atoi(ph.value);	
	float a = atof(pa.value);
	

    //! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	//! Process
	libUSTG::cflimage output = input.padding(w, h, a);

	
	//! Saving output
	output.save(pout.value); 
	
    
    return EXIT_SUCCESS;
}
