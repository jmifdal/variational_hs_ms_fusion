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
	ParStruct pinput = {"image", NULL, "input image1"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output image2"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_histogram_midway_sequence", "histogram midway equalization",
                      argc, argv, options, parameters))
		return 0;
	
	
	
	
	
	// input
    libUSTG::cflmovie input;
	input.load(pinput.value);
	
	int nframes = input.n();
	

	
	// memory
    libUSTG::cflimage *idata = new libUSTG::cflimage[nframes];
    libUSTG::cflimage *odata = new libUSTG::cflimage[nframes];
	
	
	
	for (int ii=0; ii < nframes ; ii++)
	{
		idata[ii] = input.getframe(ii);
		odata[ii] = idata[ii];
				
	}
	
	
	int width = idata[0].w();
	int height = idata[0].h();
	
	
	// process for each channel
	printf("input movie  n: %d  w: %d  h: %d \n", nframes, width, height);
	printf("dealing with channel:  ");
	for (int d = 0; d < idata[0].c(); d++)
	{
		printf("%d  ", d);
		fflush(stdout);
		
		
		// take pointers
		float ** iptr = new float*[nframes];
		float ** optr = new float*[nframes];
	
		
		for (int ii=0; ii < nframes; ii++)
		{
			iptr[ii] = idata[ii].v(d);
			optr[ii] = odata[ii].v(d);
			
		}
		

		// process
		libUSTG::fk_histogram_midway_sequence(iptr, optr, nframes, width, height);
			
		
	}
	
	printf("\n");
		
    
	//! save movie
	libUSTG::cflmovie output(pout.value, nframes);
	for (int ii=0; ii < nframes; ii++)
	{
		output.write(odata[ii]);
	}
	

}