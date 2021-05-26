#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;



int main(int argc, char **argv)
{
    vector <OptStruct *> options;
    OptStruct oPatchSize = {"p:", 0, "3", NULL, "half-size of patches for local equalization"}; options.push_back(&oPatchSize);
    OptStruct oStep = {"s:", 0, "1", NULL, "step between patches (at most, 2*patchSize+1"}; options.push_back(&oStep);
	
	vector<ParStruct *> parameters;
	ParStruct pref = {"reference", NULL, "reference input image"}; parameters.push_back(&pref);
	ParStruct pin = {"modified", NULL, "input image to be modified"}; parameters.push_back(&pin);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	
	if(!parsecmdline("src_histogram_specification_local", "Local histogram specification", argc, argv, options, parameters))
		return EXIT_FAILURE;

    
    // Parameters
    int patchSize = atoi(oPatchSize.value);
    int step = atoi(oStep.value);
    
    
    // Reference and modified images
    libUSTG::cflimage ref, inp;
    
    ref.load(pref.value);
    inp.load(pin.value);
    
    int width = ref.w();
    int height = ref.h();
    int dim = width * height;
    int num_channels = ref.c();
    
    if((inp.w() != width) || (inp.h() != height))
    {
        printf("ERROR :: src_histogram_specification_local :: reference and modified images sizes don't match.\n");
        return EXIT_FAILURE;
    }
    
    if(inp.c() != num_channels)
    {
        printf("ERROR :: src_histogram_specification_local :: reference and modified images must have the same number of channels.\n");
        return EXIT_FAILURE;
    }
    
    float **reference = new float*[num_channels];
    float **input = new float*[num_channels];
    
    for(int c = 0; c < num_channels; c++)
    {
        reference[c] = new float[dim];
        input[c] = new float[dim];
        
        libUSTG::fpCopy(ref.v(c), reference[c], dim);
        libUSTG::fpCopy(inp.v(c), input[c], dim);
    }
    
    
    // Apply local histogram specification
    if(step > 2 * patchSize + 1)
    {
        printf("ERROR :: src_histogram_specification_local :: step between patches must be at most 2*patchSize+1.\n");
        return EXIT_FAILURE;
    }
    
    int patchDim = (2 * patchSize + 1) * (2 * patchSize + 1);
    
    float *mask = new float[dim];
    float **output = new float*[num_channels];
    
    for(int c = 0; c < num_channels; c++)
    {
        output[c] = new float[dim];
        
        libUSTG::fpClear(output[c], 0.0f, dim);
        libUSTG::fpClear(mask, 0.0f, dim);
        
        for(int y = patchSize; y < height; y += step)
        {
            for(int x = patchSize; x < width; x += step)
            {
                // Adapt patch size to image boundaries
                int jmin = MAX(y - patchSize, 0);
                int jmax = MIN(y + patchSize, height - 1);
                int imin = MAX(x - patchSize, 0);
                int imax = MIN(x + patchSize, width - 1);
                
                int patchWidth = imax - imin + 1;
                int patchHeight = jmax - jmin + 1;
                int patchDim = patchWidth * patchHeight;
                
                float *patchRef = new float[patchDim];
                float *patchInp = new float[patchDim];
                float *patchOut = new float[patchDim];
                
                for(int j = jmin, jp = 0; j <= jmax; j++, jp++)
                {
                    for(int i = imin, ip = 0; i <= imax; i++, ip++)
                    {
                        int p = j * width + i;
                        int pp = jp * patchWidth + ip;
                        
                        patchRef[pp] = reference[c][p];
                        patchInp[pp] = input[c][p];
                    }
                }
                
                libUSTG::fk_histogram_specification(patchRef, patchInp, patchOut, patchWidth, patchHeight, patchWidth, patchHeight);
                
                for(int j = jmin, jp = 0; j <= jmax; j++, jp++)
                {
                    for(int i = imin, ip = 0; i <= imax; i++, ip++)
                    {
                        int p = j * width + i;
                        int pp = jp * patchWidth + ip;
                        
                        output[c][p] += patchOut[pp];
                        mask[p] = mask[p] + 1.0f;
                    }
                }
                
                
                delete[] patchRef;
                delete[] patchInp;
                delete[] patchOut;

                
            }
        }
        
        // Patch overlapping
        for(int p = 0; p < dim; p++)
            if(mask[p] != 0)
                output[c][p] /= mask[p];
    }
    
    
    // Saving output
    libUSTG::cflimage out(width, height, num_channels);
    
    for(int c = 0; c < num_channels; c++)
        libUSTG::fpCopy(output[c], out.v(c), dim);
    
    out.save(pout.value);
    
    
    // Delete allocated memory
    for(int c = 0; c < num_channels; c++)
    {
        delete[] reference[c];
        delete[] input[c];
        delete[] output[c];
    }
    
    delete[] reference;
    delete[] input;
    delete[] output;
    delete[] mask;

	
    return EXIT_SUCCESS;
}
