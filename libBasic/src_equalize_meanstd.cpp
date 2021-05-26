#include "../library/libImage.h"
#include "../library/libBasic.h"



int main(int argc, char **argv)
{
    
    std::vector <OptStruct *> options;
	
	std::vector<ParStruct *> parameters;
    ParStruct pref = {"reference", NULL, "reference input image"}; parameters.push_back(&pref);
    ParStruct pin = {"modified", NULL, "input image to be modified"}; parameters.push_back(&pin);
	ParStruct pout = {"output", NULL, "equalized image"}; parameters.push_back(&pout);
	
    if(!parsecmdline("src_equalize_meanstd", "Global mean and standard deviation equalization applied to each channel",
                     argc, argv, options, parameters))
        return EXIT_FAILURE;

    
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
        printf("ERROR :: src_equalize_meanstd :: reference and modified images sizes must match.\n");
        return EXIT_FAILURE;
    }
    
    if(inp.c() != num_channels)
    {
        printf("ERROR :: src_equalize_meanstd :: reference and modified images must have the same number of channels.\n");
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
    
    
    // Apply global equalization
    float **output = new float*[num_channels];
    
    for(int c = 0; c < num_channels; c++)
    {
        output[c] = new float[dim];
        
        float mean_inp = 0.0f;
        float std_inp = 0.0f;
        float mean_ref = 0.0f;
        float std_ref = 0.0f;
        
        for(int i = 0; i < dim; i++)
        {
            mean_inp += input[c][i];
            mean_ref += reference[c][i];
            
            std_inp += (input[c][i] * input[c][i]);
            std_ref += (reference[c][i] * reference[c][i]);
        }
        
        float fdim = (float) dim;
        
        mean_inp /= fdim;
        mean_ref /= fdim;
        std_inp /= fdim;
        std_ref /= fdim;
        
        std_inp = sqrtf(fabs(std_inp - mean_inp * mean_inp));
        std_ref = sqrtf(fabs(std_ref - mean_ref * mean_ref));
        
        float std_div = std_ref / std_inp;
        
        for(int i = 0; i < dim; i++)
            output[c][i] = mean_ref + std_div * (input[c][i] - mean_inp);
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
    

	return EXIT_SUCCESS;
}
