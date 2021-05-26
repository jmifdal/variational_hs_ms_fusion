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
    ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
    ParStruct poutput = {"output", NULL, "output image"}; parameters.push_back(&poutput);
    ParStruct pV = {"V", NULL, "save right-singular matrix"}; parameters.push_back(&pV);
    
    if(!parsecmdline("src_pca_color_transform_back", "Apply PCA back to original color space", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Input image
    libUSTG::cflimage in;
    in.load(pinput.value);
    
    int width = in.w();
    int height = in.h();
    int dim = width * height;
    int num_channels = in.c();
    
    float **input = new float*[num_channels];
    for(int c = 0; c < num_channels; c++)
    {
        input[c] = new float[dim];
        libUSTG::fpCopy(in.v(c), input[c], dim);
    }
    
    
    // Input V^T vector
    libUSTG::flimage Vin;
    Vin.load(pV.value);
    
    if((Vin.w() != num_channels) || (Vin.h() != num_channels))
    {
        printf("ERROR :: src_pca_color_transform_back :: dimensions of V must coincide with number of image channels.\n");
        return EXIT_FAILURE;
    }
    
    float **V = new float*[num_channels];
    
    for(int c1 = 0; c1 < num_channels; c1++)
    {
        V[c1] = new float[num_channels];
        
        for(int c2 = 0; c2 < num_channels; c2++)
            V[c1][c2] = Vin[c2 * num_channels + c1];
    }
    
    
    // Apply PCA transformation back. Force the first component to be negative.
    float **output = new float*[num_channels];
    for(int c = 0; c < num_channels; c++)
        output[c] = new float[dim];
    
    for(int i = 0; i < dim; i++)
    {
        for(int c1 = 0; c1 < num_channels; c1++)
        {
            output[c1][i] = 0.0f;
            
            for(int c2 = 0; c2 < num_channels; c2++)
            {
                if(c2 == 0)
                    output[c1][i] += (-input[c2][i] * V[c2][c1]);
                else
                    output[c1][i] += (input[c2][i] * V[c2][c1]);
            }
        }
    }

    
    // Save output
    libUSTG::cflimage out(width, height, num_channels);
    
    for(int c = 0; c < num_channels; c++)
        libUSTG::fpCopy(output[c], out.v(c), dim);
        
    out.save(poutput.value);
    
    
    // Delete allocated memory
    for(int c = 0; c < num_channels; c++)
    {
        delete[] input[c];
        delete[] output[c];
        delete[] V[c];
    }
    
    delete[] input;
    delete[] output;
    delete[] V;
    
    
    return EXIT_SUCCESS;
}


