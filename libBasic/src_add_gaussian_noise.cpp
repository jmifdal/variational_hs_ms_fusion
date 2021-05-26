#include "../library/libBasic.h"
#include "../library/libImage.h"


int main(int argc, char **argv)
{
	
    std::vector <OptStruct *> options;
	OptStruct ostd = {"s:", 0, NULL, NULL, "noise standard deviation"}; options.push_back(&ostd);
    OptStruct osnr = {"r:", 0, NULL, NULL, "signal to noise ratio"}; options.push_back(&osnr);
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"output", NULL, "output noisy image"}; parameters.push_back(&pout);
	
    char *arg = "INFO: Add gaussian noise of given standard deviation or estimate it according to given signal to noise ratio";
    
    if(!parsecmdline("src_add_gaussian_noise", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    
    // Process
    if(ostd.flag)
    {
        // Noise standard deviation
        float std = atof(ostd.value);
        
        // Add Gaussian noise
        input.addGaussianNoise(std);
        
        // Save output
        input.save(pout.value);
        
    } else
    {
        if(osnr.flag)
        {
            // Input image
            int width = input.w();
            int height = input.h();
            int num_channels = input.c();
            int dim = width * height;
            
            float **image = new float*[num_channels];
            float **noisy = new float*[num_channels];
            
            for(int c = 0; c < num_channels; c++)
            {
                image[c] = new float[dim];
                libUSTG::fpCopy(input.v(c), image[c], dim);
                
                noisy[c] = new float[dim];
            }
            
            // Compute noise standard deviation according to SNR of each spectral band
            float snr = atof(osnr.value);
            float *sd_values = new float[num_channels];
            libUSTG::compute_noise_sd(image, sd_values, width, height, num_channels, snr);
            
            // Add Gaussian noise
            for(int c = 0; c < num_channels; c++)
                libUSTG::fpAddNoiseGaussian(image[c], noisy[c], sd_values[c], 0, dim);
            
            // Save output
            libUSTG::cflimage output(width, height, num_channels);
            
            for(int c = 0; c < num_channels; c++)
                libUSTG::fpCopy(noisy[c], output.v(c), dim);
            
            output.save(pout.value);
            
            // Free memory
            for(int c = 0; c < num_channels; c++)
            {
                delete[] image[c];
                delete[] noisy[c];
            }
            
            delete[] image;
            delete[] noisy;
            delete[] sd_values;
            
        } else
        {
            printf("ERROR :: src_add_gaussian_noise :: either noise standard deviation or signal to noise ratio must be provided.\n");
            return EXIT_FAILURE;
        }
    }
    
    
    return EXIT_SUCCESS;
}
