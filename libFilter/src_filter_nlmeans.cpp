#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"
#include "../library/libFilteringNLmeans.h"

int main(int argc, char **argv)
{
	
    std::vector <OptStruct *> options;
	OptStruct obloc = {"b:", 0, "15", NULL, "research block size"};  options.push_back(&obloc);	
	OptStruct owin  = {"w:", 0, "3", NULL, "comparison window size"}; options.push_back(&owin);	
	OptStruct otype  = {"t:", 0, "0", NULL, "type"}; options.push_back(&otype);
	OptStruct osecimage = 	 {"i:", 0, NULL, NULL, "optional second image for window distance computation"}; options.push_back(&osecimage);
	OptStruct omulti = 	{"h:", 0, "1.0", NULL, "multiplier for image filtering"}; options.push_back(&omulti);
	OptStruct ospatial = {"s:", 0, "0.0", NULL, "multiplier for spatial filtering"}; options.push_back(&ospatial);
	OptStruct omask = {"m:", 0, NULL, NULL, "initial mask"}; options.push_back(&omask);
	OptStruct oIter = {"n:", 0, "1", NULL, "number of iterations"}; options.push_back(&oIter);
	
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_filter_nlmeans", "Filters by nonlocal means", argc, argv, options, parameters))
		return 0;
	
	
	int bloc = (atoi(obloc.value) - 1) / 2;
	int win = (atoi(owin.value) - 1) / 2;
	
	float  fImFilt =  atof(omulti.value);
	float  fSpaFilt =  atof(ospatial.value);
	int nIter = atoi(oIter.value);
    
	
	
	//! Input & output
    libUSTG::cflimage input;
	input.load(pinput.value);
	
	int width = input.w();
	int height = input.h();
	int nchannels = input.c();	
	
	
	float **fpI = new float*[input.c()];
	float **fpO = new float*[input.c()];
	float **fpA = NULL;
	
    libUSTG::cflimage output(width, height, nchannels);
	for (int ii=0; ii < input.c(); ii++) {
		
		fpI[ii] = input.v(ii);
		fpO[ii] = output.v(ii);
		
	}
	
	
	//! Secondary image
	libUSTG::cflimage input2;
    int iAchannels = 0;
	if (osecimage.flag)
	{
		
		input2.load(osecimage.value);
		if (input.w() != input2.w() || input.h() != input2.h()) 
		{		
			printf("secondary image of different size\n");
			exit(-1);
		}
		
        iAchannels = input2.c();
		fpA = new float*[iAchannels];
        for (int ii=0; ii < iAchannels; ii++) fpA[ii] = input2.v(ii);
	}
	
    
    
    //! Mask image
    libUSTG::cflimage input3;
    float *fpM = NULL;
    if (omask.flag)
    {
        input3.load(omask.value);
        fpM = input3.v();
        
    } else
    {
        fpM = new float[width*height];
        libUSTG::fpClear(fpM, 1.0, width*height);
    }
    
    
    if (nIter>1)
        libUSTG::fpNeighFilteringIt(
                    win,                // Half size of comparison window
                    bloc,               // Half size of research window
                    fSpaFilt,            // Filtering parameter
                    fImFilt,            // Filtering parameter
                    nIter,
                    fpI,                // Input
                    fpA,                // Auxiliary input for distance computation
                    fpO,                // Output
                    fpM,
                    atoi(otype.value),
                    nchannels, iAchannels,  width, height);
    else
        libUSTG::fpNeighFiltering(
                                  win,                // Half size of comparison window
                                  bloc,               // Half size of research window
                                  fSpaFilt,            // Filtering parameter
                                  fImFilt,            // Filtering parameter
                                  fpI,                // Input
                                  fpA,                // Auxiliary input for distance computation
                                  fpO,                // Output
                                  fpM,
                                  atoi(otype.value),
                                  nchannels, iAchannels,  width, height);
    
        
        
    if (!omask.flag) delete[] fpM;
    if (fpA != NULL) delete[] fpA;
    delete[] fpI;
    delete[] fpO;
    
	output.save(pout.value);
    
}



