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
    ParStruct pInput = {"input", NULL, "input ground-truth high-resolution hyperspectral image"}; parameters.push_back(&pInput);
    ParStruct pOutput = {"output", NULL, "output high-resolution multispectral image"}; parameters.push_back(&pOutput);
    ParStruct pOperator = {"operator", NULL, "spectral-downsampling matrix operator"}; parameters.push_back(&pOperator);
   
    char* arg = "INFO: Generate high-resolution multispectral image from ground-truth high-resolution hyperspectral image";
    
    if(!parsecmdline("src_hyperspectral_generate_multispectral", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Ground-truth high-resolution hyperspectral image
    libUSTG::cflimage input;
    input.load(pInput.value);
    
    int width = input.w();
    int height = input.h();
    int dim = width * height;
    int hs_channels = input.c();
    
    float **truth = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        truth[h] = new float[dim];
        libUSTG::fpCopy(input.v(h), truth[h], dim);
    }
    
    
    // Spectral-downsampling operator
    libUSTG::flimage Sinput;
    Sinput.load(pOperator.value);
    
    if(Sinput.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_generate_multispectral :: spectral-downsampling matrix must be a one-channel image\n");
        return EXIT_FAILURE;
    }
    
    if(Sinput.w() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_generate_multispectral :: spectral-downsampling matrix width does not coincide with number of hyperspectral channels\n");
        return EXIT_FAILURE;
    }
    
    int ms_channels = Sinput.h();
    
    
    float **Smatrix = new float*[ms_channels];
    
    for(int m = 0; m < ms_channels; m++)
    {
        Smatrix[m] = new float[hs_channels];
        
        for(int h = 0; h < hs_channels; h++)
            Smatrix[m][h] = Sinput.v()[m * hs_channels + h];
    }
    
    
    // Compute high-resolution multispectral image
    float **multispectral = new float*[ms_channels];
    for(int m = 0; m < ms_channels; m++)
        multispectral[m] = new float[dim];
    
    for(int m = 0; m < ms_channels; m++)
    {
        for(int i = 0; i < dim; i++)
        {
            float Sm = 0.0f;
            float Sw = 0.0f;
            
            for(int h = 0; h < hs_channels; h++)
            {
                Sm += Smatrix[m][h] * truth[h][i];
                Sw += Smatrix[m][h];
            }
            
            multispectral[m][i] = Sm / Sw;
        }
    }
    
    
    // Save output
    libUSTG::cflimage output(width, height, ms_channels);
    
    for(int m = 0; m < ms_channels; m++)
        libUSTG::fpCopy(multispectral[m], output.v(m), dim);
    
    output.save(pOutput.value);
    
    
    // Delete allocated memory
    for(int h = 0; h < hs_channels; h++)
        delete[] truth[h];
    
    for(int m = 0; m < ms_channels; m++)
    {
        delete[] Smatrix[m];
        delete[] multispectral[m];
    }
    
    delete[] truth;
    delete[] Smatrix;
    delete[] multispectral;
    
    return EXIT_SUCCESS;
}
