#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "../library/libImage.h"


int main(int argc, char **argv)
{
    
    std::vector <OptStruct *> options;
    
    std::vector<ParStruct *> parameters;
    ParStruct pInput = {"input", NULL, "input  multispectral image"}; parameters.push_back(&pInput);
    ParStruct pOutput = {"output", NULL, "output hyperspectral-pan image"}; parameters.push_back(&pOutput);
    ParStruct pOperator = {"operator", NULL, "spectral-downsampling matrix operator"}; parameters.push_back(&pOperator);
   
    char* arg = "INFO: Generate hyperspectral pan image as linear combination of multispectral channels";
    
    if(!parsecmdline("src_hyperspectral_generate_hyperpan", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Multispectral image
    libUSTG::cflimage input;
    input.load(pInput.value);
    
    int width = input.w();
    int height = input.h();
    int dim = width * height;
    int ms_channels = input.c();
    
    float **multispectral = new float*[ms_channels];
    
    for(int m = 0; m < ms_channels; m++)
    {
        multispectral[m] = new float[dim];
        libUSTG::fpCopy(input.v(m), multispectral[m], dim);
    }
    
    
    // Spectral-downsampling operator
    libUSTG::flimage Sinput;
    Sinput.load(pOperator.value);
    
    if(Sinput.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_generate_hyperpan :: spectral-downsampling matrix must be a one-channel image\n");
        return EXIT_FAILURE;
    }
    
    if(Sinput.h() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_generate_hyperpan :: spectral-downsampling matrix height does not coincide with number of multispectral channels\n");
        return EXIT_FAILURE;
    }
    
    int hs_channels = Sinput.w();
    
    float **Smatrix = new float*[ms_channels];
    
    for(int m = 0; m < ms_channels; m++)
    {
        Smatrix[m] = new float[hs_channels];
        
        for(int h = 0; h < hs_channels; h++)
            Smatrix[m][h] = Sinput.v()[m * hs_channels + h];
    }
    
    
    // Compute hyperspectral-pan image
    float **hyper_pan = new float*[hs_channels];
    
    for(int h = 0; h < hs_channels; h++)
        hyper_pan[h] = new float[dim];
    
    for(int h = 0; h < hs_channels; h++)
    {
        for(int i = 0; i < dim; i++)
        {
            float Sm = 0.0f;
            float Sw = 0.0f;
            
            for(int m = 0; m < ms_channels; m++)
            {
                Sw += Smatrix[m][h];
                Sm += Smatrix[m][h] * multispectral[m][i];
            }
            
            hyper_pan[h][i] = Sm / Sw;
        }
    }
    
    
    // Save output
    libUSTG::cflimage output(width, height, hs_channels);
    
    for(int h = 0; h < hs_channels; h++)
        libUSTG::fpCopy(hyper_pan[h], output.v(h), dim);
    
    output.save(pOutput.value);
    
    
    // Delete allocated memory
    for(int h = 0; h < hs_channels; h++)
        delete[] hyper_pan[h];
    
    for(int m = 0; m < ms_channels; m++)
    {
        delete[] Smatrix[m];
        delete[] multispectral[m];
    }
    
    delete[] multispectral;
    delete[] Smatrix;
    delete[] hyper_pan;
    
    return EXIT_SUCCESS;
}
