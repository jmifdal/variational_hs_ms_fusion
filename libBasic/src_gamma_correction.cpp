




#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;


int main(int argc, char **argv)
{
	

	
	vector <OptStruct *> options;
    OptStruct oM = {"m:", 0,  "255", NULL,"max value"};  options.push_back(&oM);

	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct poutput = {"image", NULL, "output image"}; parameters.push_back(&poutput);
    ParStruct pgamma = {"gamma", NULL, "gamma value"}; parameters.push_back(&pgamma);
	
	
	if (!parsecmdline("src_gamma_correction", "applies gamma correction", argc, argv, options, parameters))
		return 0;
	


	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
    
	
	//! Process
    float fMax = atof(oM.value);
    float gamma = atof(pgamma.value);
    for (int ii=0; ii < input.whc(); ii++)
    {
        input[ii] = fMax * powf(input[ii] / fMax, gamma);
    }
	
    
    //! Save
	input.save(poutput.value);
	return 1;
	
}
