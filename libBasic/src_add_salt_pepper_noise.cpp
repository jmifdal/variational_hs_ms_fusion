#include "../library/libBasic.h"
#include "../library/libImage.h"


int main(int argc, char **argv)
{
    
    std::vector <OptStruct *> options;
    OptStruct op = {"p:", 0, "0.15", NULL,"amount of noisy pixels"}; options.push_back(&op);
    
    std::vector<ParStruct *> parameters;
    ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
    ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
    
    
    if(!parsecmdline("src_add_salt_pepper_noise","Adds salt-and-pepper noise to (p*10)% of pixels",
                      argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    // Parameters
    float p = atof(op.value);
    float p2 = p / 2.0f;
    float p3 = 1.0f - p2;
    
    // Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    // Image sizes
    int width = input.w();
    int height = input.h();
    int dim = width * height;
    int num_channels = input.c();
    
    // Add salt-and-pepper noise
    float **output = new float*[num_channels];
    
    for(int c = 0; c < num_channels; c++)
    {
        output[c] = new float[dim];
        
        for(int i = 0; i < dim; i++)
        {
            float tmp = (float) rand() / RAND_MAX;
            
            if(tmp <= p2)
                output[c][i] = 0.0f;
            else if(tmp >= p3)
                output[c][i] = 255.0f;
            else
                output[c][i] = input.v(c)[i];
        }
    }
    
    // Save result
    libUSTG::cflimage out(width, height, num_channels);
    
    for(int c = 0; c < num_channels; c++)
        libUSTG::fpCopy(output[c], out.v(c), dim);
    
    out.save(pout.value);
    
    // Delete allocated memory
    for(int c = 0; c < num_channels; c++)
        delete[] output[c];
    
    delete[] output;
    
    return EXIT_SUCCESS;
}
