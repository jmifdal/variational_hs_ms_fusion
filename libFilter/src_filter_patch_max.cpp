#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"
#include "../library/libFilteringNLmeans.h"

int main(int argc, char **argv)
{
	
    std::vector <OptStruct *> options;
	OptStruct oR  =   {"r:", 0, "2.0", NULL,"radius"}; options.push_back(&oR);
	OptStruct omask = 	{"m:", 0, NULL, NULL, "initial mask"}; options.push_back(&omask);
	OptStruct oIter = 	{"n:", 0, "1", NULL, "number of iterations"}; options.push_back(&oIter);
	
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_filter_patch_max", "Filters patches by maximum", argc, argv, options, parameters))
		return 0;
	
	

	//! Input & output
    libUSTG::cflimage input;
	input.load(pinput.value);
	
	
    libUSTG::cflimage output(input.w(), input.h(), input.c());
    output = input;
    
    float *imMask = NULL;
    if (omask.flag)
    {
        libUSTG::flimage tmpMask; tmpMask.load(omask.value);

        imMask = new float[tmpMask.wh()];
        libUSTG::fpCopy(tmpMask.v(), imMask, tmpMask.wh());
        
    }
    
    
    int nIter = atoi(oIter.value);
    float fRadius=atof(oR.value);
    
    
    for (int i=0; i < nIter; i++)
    {
    
        
        for (int iC=0; iC < input.c(); iC++)
        {
            libUSTG::fiPatchMax(input.v(iC), output.v(iC), fRadius, input.w(), input.h());
            
        }
 
        if (imMask)
        {
            for (int ii=0; ii < input.wh(); ii++)
            {
                if (imMask[ii] > 0.0f)
                    for (int iC=0; iC < input.c(); iC++)
                    {
                        input[iC * input.c() + ii] = output[iC * input.c() + ii] ;
                    }
                
            }
            
        }
        else
            input=output;
    }
    
    
    if (imMask) delete[] imMask;
	output.save(pout.value);
    
    
}



