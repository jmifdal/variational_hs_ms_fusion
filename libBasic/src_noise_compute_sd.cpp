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
    ParStruct pout = {"sd", NULL, "sd values"}; parameters.push_back(&pout);
    ParStruct pSNR = {"SNR", NULL, "signal to noise ratio"}; parameters.push_back(&pSNR);
    
    if(!parsecmdline("src_noise_compute_sd", "Compute the noise standard deviation of each spectral band according to SNR",
                     argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Parameters
    float snr = atof(pSNR.value);
    
    
    // Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    int width = input.w();
    int height = input.h();
    int num_channels = input.c();
    int dim =  width * height;
    
    float **image = new float*[num_channels];
    
    for(int c = 0; c < num_channels; c++)
    {
        image[c] = new float[dim];
        libUSTG::fpCopy(input.v(c), image[c], dim);
    }
    
    
    // Process
    float *sd_values = new float[num_channels];
    libUSTG::compute_noise_sd(image, sd_values, width, height, num_channels, snr);
    
    
    // Output
    libUSTG::flimage output(num_channels, 1);
    libUSTG::fpCopy(sd_values, output.v(), dim);
    output.save(pout.value);
    
    
    // Free memory
    for(int c = 0; c <num_channels; c++)
        delete[] image[c];
    
    delete[] image;
    delete[] sd_values;
    
    
    return EXIT_SUCCESS;
}

