#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "../library/libImage.h"
#include "../library/libFilteringVariational.h"


int main(int argc, char **argv)
{
    
    std::vector <OptStruct *> options;
    OptStruct oLmbTV = {"l:", 0, "1", NULL, "Trade-off parameter for TV regularization"}; options.push_back(&oLmbTV);
    OptStruct oLmbNL = {"b:", 0, "1", NULL, "Trade-off parameter for NL regularization"}; options.push_back(&oLmbNL);
    OptStruct ohPatch = {"p:", 0, "2.5", NULL, "Filtering parameter for patch-similarity distance in NL weights"}; options.push_back(&ohPatch);
    OptStruct ohSpatial = {"s:", 0, "2.5", NULL, "Filtering parameter for spatial distance in NL weights"}; options.push_back(&ohSpatial);
    OptStruct onNeigh = {"n:", 0, "25", NULL, "Number of neighbours used in NL weights"}; options.push_back(&onNeigh);
    OptStruct oreswind = {"r:", 0, "10", NULL, "Half-size of research window"}; options.push_back(&oreswind);
    OptStruct ocompwind = {"c:", 0, "1", NULL, "Half-size of comparison window"}; options.push_back(&ocompwind);
    OptStruct oflagNormW = {"w", 0, NULL, NULL, "Normalize weights"};  options.push_back(&oflagNormW);
    OptStruct oFlagUseMaskW = {"u", 0, NULL, NULL, "Use binary mask to compute weights"}; options.push_back(&oFlagUseMaskW);
    OptStruct oTol = {"e:", 0, "0.000001", NULL, "Stopping precision"}; options.push_back(&oTol);
    OptStruct oMaxIter = {"m:", 0, "1000", NULL, "max algorithm iterations"}; options.push_back(&oMaxIter);
    OptStruct oTau = {"t:", 0, "0.05", NULL, "step-size of proximity operator of data-fidelity term)"}; options.push_back(&oTau);
    OptStruct oSigma = {"g:", 0, "0.05", NULL, "step-size of proximity operator of regularization term)"}; options.push_back(&oSigma);
    
    std::vector<ParStruct *> parameters;
    ParStruct pInput = {"input", NULL, "initialization"}; parameters.push_back(&pInput);
    ParStruct pOutput = {"output", NULL, "filtered output"}; parameters.push_back(&pOutput);
    ParStruct pWimage = {"weight_image", NULL, "image on which nonlocal weights are computed"}; parameters.push_back(&pWimage);
    ParStruct pMask = {"mask", NULL, "binary mask for activated pixels in the data-fidelity term"}; parameters.push_back(&pMask);
    
    char* arg = "INFO: Filtering an image by the variational model composed of TV + NL regularizations and weighted L2 fidelity term.\n"
                "      The NL-regularization weights are computed on the provided image using bilateral weights (patch-similarity + spatial distance).\n"
                "      The L2 fidelity term to the initialization is weighted by the provided binary mask.";
    
    if(!parsecmdline("src_filter_variational_TVNL_L2w", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Parameters
    float lmbTV = atof(oLmbTV.value);
    float lmbNL = atof(oLmbNL.value);
    float hPatch = atof(ohPatch.value);
    float hSpatial = atof(ohSpatial.value);
    int numNeigh = atoi(onNeigh.value);
    int reswind = atoi(oreswind.value);
    int compwind = atoi(ocompwind.value);
    int flagNormalizeW = oflagNormW.flag;
    int flagUseMaskW = oFlagUseMaskW.flag;
    float tol = atof(oTol.value);
    int maxIter = atoi(oMaxIter.value);
    float tau = atof(oTau.value);
    float sigma = atof(oSigma.value);
    
    
    // Initialization
    libUSTG::flimage input;
    input.load(pInput.value);
    
    int width = input.w();
    int height = input.h();
    int dim = width * height;
    
    float *f = new float[dim];
    libUSTG::fpCopy(input.v(), f, dim);
    
    
    // Image on which nonlocal weights are computed
    libUSTG::cflimage img;
    img.load(pWimage.value);
    
    if((img.w() != width) || (img.h() != height))
    {
        printf("ERROR :: src_filter_variational_TVNL_L2w :: dimensions of weight_image and initialization don't match\n");
        return EXIT_FAILURE;
    }
    
    int num_channels = img.c();
    
    float **wimage = new float*[num_channels];
    
    for(int c= 0; c < num_channels; c++)
    {
        wimage[c] = new float[dim];
        libUSTG::fpCopy(img.v(c), wimage[c], dim);
    }
    
    
    // Binary mask for weighted L2 data term
    libUSTG::flimage fmask;
    fmask.load(pMask.value);
    
    if((fmask.w() != width) || (fmask.h() != height))
    {
        printf("ERROR :: src_filter_variational_TVNL_L2w :: dimensions of mask and initialization don't match\n");
        return EXIT_FAILURE;
    }
    
    float *mask = new float[dim];
    libUSTG::fpCopy(fmask.v(), mask, dim);
    
    
    // Compute NLTV weights
    std::vector< std::vector<float> > wxy;
    std::vector< std::vector<int> > posxy;
    std::vector< std::vector<float> > wyx;
    std::vector< std::vector<int> > posyx;
    std::vector< std::vector<int> > posw;
    
    if(flagUseMaskW)
    {
        libUSTGFiltVar::nlweights_reg_mask(wimage, mask, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeigh,
                                           reswind, compwind, flagNormalizeW, num_channels, width, height);
        
    } else {
        
        libUSTGFiltVar::nlweights_reg(wimage, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeigh,
                                      reswind, compwind, flagNormalizeW, num_channels, width, height);
    }
    

    // Variational filtering
    float *u = new float[dim];
    
    libUSTGFiltVar::filtering_TV_NL_L2w(u, f, mask, wxy, posxy, wyx, posyx, posw, lmbTV,
                                        lmbNL, tau, sigma, tol, maxIter, width, height);
    
    
    // Save output
    libUSTG::flimage output(width, height, u);
    output.save(pOutput.value);
    
    
    // Delete allocated memory
    for(int c = 0; c < num_channels; c++)
        delete[] wimage[c];
    
    delete[] f;
    delete[] wimage;
    delete[] mask;
    delete[] u;
    
    return EXIT_SUCCESS;
}
