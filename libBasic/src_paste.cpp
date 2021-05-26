#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input1", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"input2", NULL, "input image"}; parameters.push_back(&pinput2);
	ParStruct pout = {"output", NULL, "output image"}; parameters.push_back(&pout);
	ParStruct px = {"x", NULL, "top left x coordinate"}; parameters.push_back(&px);
	ParStruct py = {"y", NULL, "top left y coordinate"}; parameters.push_back(&py);
	 
	if(!parsecmdline("src_paste", "Paste image inside another one with left top corner (x,y)",
                     argc, argv, options, parameters))
        return EXIT_FAILURE;
	

	//! Parameters
	int x = atoi(px.value);
	int y = atoi(py.value);
	
	
	//! Inputs
	libUSTG::cflimage input, input2;
	input.load(pinput.value);
	input2.load(pinput2.value);
	
    
	//! Process
	input.paste(input2, x, y);
	
    
	//! Saving output
	input.save(pout.value);
    
    
	return EXIT_SUCCESS;
}
