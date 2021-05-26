#include "../library/libImage.h"
#include "../library/libBasic.h"



int main(int argc, char **argv)
{

    std::vector <OptStruct *> options;
    OptStruct oPatchSize = {"p:", 0, "3", NULL, "half-size of patches for local equalization"}; options.push_back(&oPatchSize);
    OptStruct oStep = {"s:", 0, "1", NULL, "step between patches (at most, 2*patchSize+1"}; options.push_back(&oStep);
	
    std::vector<ParStruct *> parameters;
    ParStruct pref = {"reference", NULL, "reference input image"}; parameters.push_back(&pref);
    ParStruct pin = {"modified", NULL, "input image to be modified"}; parameters.push_back(&pin);
	ParStruct pout = {"output", NULL, "equalized image"}; parameters.push_back(&pout);
	
    if(!parsecmdline("src_equalize_meanstd_local", "Local mean and standard deviation equalization applied to each channel",
                     argc, argv, options, parameters))
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
        printf("ERROR :: src_equalize_meanstd_local :: input and reference images sizes don't match.\n");
        return EXIT_FAILURE;
    }
    
    if(inp.c() != num_channels)
    {
        printf("ERROR :: src_equalize_meanstd_local :: input and reference images must have the same number of channels.\n");
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
    
    
    // Apply local equalization
    if(step > 2 * patchSize + 1)
    {
        printf("ERROR :: src_equalize_meanstd_local :: step between patches must be at most 2*patchSize+1.\n");
        return EXIT_FAILURE;
    }
    
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
                
                // Compute mean and standard deviations on current patch
                float mean_inp = 0.0f;
                float std_inp = 0.0f;
                float mean_ref = 0.0f;
                float std_ref = 0.0f;
                int numPixels = 0;
                
                for(int j = jmin; j <= jmax; j++)
                {
                    for(int i = imin; i <= imax; i++)
                    {
                        int p = j * width + i;
                        numPixels++;
                        
                        mean_inp += input[c][p];
                        mean_ref += reference[c][p];
                        
                        std_inp += (input[c][p] * input[c][p]);
                        std_ref += (reference[c][p] * reference[c][p]);
                    }
                }
                
                float fnumPixels = (float) numPixels;
                
                mean_inp /= fnumPixels;
                mean_ref /= fnumPixels;
                std_inp /= fnumPixels;
                std_ref /= fnumPixels;
                
                std_inp = sqrtf(fabs(std_inp - mean_inp * mean_inp));
                std_ref = sqrtf(fabs(std_ref - mean_ref * mean_ref));
                
                float std_div = std_ref / std_inp;
                
                // Local equalization
                for(int j = jmin; j <= jmax; j++)
                {
                    for(int i = imin; i <= imax; i++)
                    {
                        int p = j * width + i;
                        output[c][p] += (mean_ref + std_div * (input[c][p] - mean_inp));
                        mask[p] = mask[p] + 1.0f;
                    }
                }
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
