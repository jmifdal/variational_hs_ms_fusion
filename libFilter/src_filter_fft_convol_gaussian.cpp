#include "../library/libBasic.h"
#include "../library/libImage.h"


int main(int argc, char **argv)
{
	
    std::vector <OptStruct *> options;
	OptStruct oStd = {"s:", 0, "1.0", NULL, "Guassian standard deviation"}; options.push_back(&oStd);
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"clear", NULL, "input clear image"}; parameters.push_back(&pinput);
	ParStruct pout = {"convolved", NULL, "output convolved image"}; parameters.push_back(&pout);
	
    
	if (!parsecmdline("src_filter_fft_convol_gaussian", "Convol image with Gaussian kernel of provided std using FFT",
                      argc, argv, options, parameters))
    {
		return EXIT_FAILURE;
    }
    
    
    //! Parameters
	float std = atof(oStd.value);
    
    
	//! Clear image
	libUSTG::cflimage input;
	input.load(pinput.value);
	
    int width = input.w();
    int height = input.h();
    int dim = width * height;
    int num_channels = input.c();
    
    
    //! Vectors
    float **data = new float*[num_channels];
    float **convol = new float*[num_channels];
    
    for(int k = 0; k < num_channels; k++)
    {
        data[k] = new float[dim];
        convol[k] = new float[dim];
        
        libUSTG::fpCopy(input.v(k), data[k], dim);
    }
	
    
    //! Convolve using FFT
    libUSTG::fft_Gaussian_convolution(convol, data, std, num_channels, width, height);
	
	
	//! Save result
    libUSTG::cflimage output(width, height, num_channels);
    
    for(int k = 0; k < num_channels; k++)
        libUSTG::fpCopy(convol[k], output.v(k), dim);
    
	output.save(pout.value);
	
    
    //! Free memory
    for(int k = 0; k < num_channels; k++)
    {
        delete[] data[k];
        delete[] convol[k];
    }
    
    delete[] data;
    delete[] convol;
    
    return EXIT_SUCCESS;
	
}