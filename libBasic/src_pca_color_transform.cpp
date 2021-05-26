#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char **argv)
{
    
    vector <OptStruct *> options;
    OptStruct oV = {"v:", 0, NULL, NULL, "save right-singular vectors"}; options.push_back(&oV);
    
    vector<ParStruct *> parameters;
    ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
    ParStruct poutput = {"output", NULL, "output image"}; parameters.push_back(&poutput);
    
    if(!parsecmdline("src_pca_color_transform", "Apply PCA to color data and provide image in decorrelated space", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Input
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
    
    
    // Compute SVD
    libUSTG::laVector S(num_channels);
    libUSTG::laMatrix US(dim, num_channels), V(num_channels, num_channels), X(dim, num_channels);
    
    for(int i = 0; i < dim; i++)
        for(int c = 0; c < num_channels; c++)
            X[i][c] = input[c][i];
    
    libUSTG::compute_pca_svd(X, S, V, US);
    
    
    // The PCA transform of input data is the image in US. Force the first component (average across channels) to be positive.
    float **output = new float*[num_channels];
    for(int c = 0; c < num_channels; c++)
        output[c] = new float[dim];
    
    for(int i = 0; i < dim; i++)
        output[0][i] = fabs(US[i][0]);
    
    for(int c = 1; c < num_channels; c++)
        for(int i = 0; i < dim; i++)
            output[c][i] = US[i][c];
    
    
    // Save right-singular vectors V^T if required
    int flagV = oV.flag;
    
    if(flagV)
    {
        libUSTG::flimage Vim(num_channels, num_channels);
        
        for(int c1 = 0; c1 < num_channels; c1++)
            for(int c2 = 0; c2 < num_channels; c2++)
                Vim[c1 * num_channels + c2] = V[c1][c2];
        
        Vim.save(oV.value);
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
    }
    
    delete[] input;
    delete[] output;
    
    
    return EXIT_SUCCESS;
}


