#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;


int main(int argc, char **argv)
{
	

	
	vector <OptStruct *> options;
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct poutput = {"image", NULL, "output image"}; parameters.push_back(&poutput);
	
	
	if (!parsecmdline("src_sqrt", "computes square root", argc, argv, options, parameters))
		return 0;
	


	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);

	
    //! Process
    for (int ii=0; ii < input.whc(); ii++) input[ii] = sqrtf(input[ii]);
    
	//! Save
	input.save(poutput.value);
	return 1;
	
}
