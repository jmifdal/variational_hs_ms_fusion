#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;


int main(int argc, char **argv)
{
    vector <OptStruct *> options;
    OptStruct oK1 = {"k:", 0, "0.01", NULL, "constant C1 = (K1*L)^2"}; options.push_back(&oK1);
    OptStruct oK2 = {"l:", 0, "0.03", NULL, "constant C2 = (K2*L)^2"}; options.push_back(&oK2);
    OptStruct oSize = {"z:", 0, "11", NULL, "size of Gaussian local window"}; options.push_back(&oSize);
    OptStruct oStd = {"s:", 0, "1.5", NULL, "Standard deviation of Gaussian local window"}; options.push_back(&oStd);
    OptStruct oRang = {"r:", 0, "255", NULL, "dynamic range of the images"};  options.push_back(&oRang);
    OptStruct oB = {"b:", 0, "0", NULL, "boundary elimination "};  options.push_back(&oB);
    
    vector<ParStruct *> parameters;
    ParStruct pinput1 = {"image1", NULL, "image1"}; parameters.push_back(&pinput1);
    ParStruct pinput2 = {"image2", NULL, "image2"}; parameters.push_back(&pinput2);
    
    if(!parsecmdline("src_ssim", "compute SSIM index", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Input images
    libUSTG::cflimage cinput1, cinput2;
    cinput1.load(pinput1.value);
    cinput2.load(pinput2.value);
    
    assert(cinput1.c() == cinput2.c() && cinput1.w() == cinput2.w() && cinput1.h() == cinput2.h());
    
    
    // Remove boundary
    int boundary = atoi(oB.value);
    
    libUSTG::cflimage input1 = cinput1.copy(boundary/2, boundary/2, cinput1.w() - boundary, cinput1.h()-boundary);
    libUSTG::cflimage input2 = cinput2.copy(boundary/2, boundary/2, cinput2.w() - boundary, cinput2.h()-boundary);
    
    cinput1.erase();
    cinput2.erase();
    
    int width = input1.w();
    int height = input1.h();
    int dim = width * height;
    int num_channels = input1.c();
    
    
    // Get grayscale images where SSIM index is computed
    float *image1 = new float[dim];
    float *image2 = new float[dim];
    
    if(num_channels == 1)
    {
        libUSTG::fpCopy(input1.v(), image1, dim);
        libUSTG::fpCopy(input2.v(), image2, dim);
        
    } else if(num_channels == 3)
    {
        for(int i = 0; i < dim; i++)
        {
            image1[i] = 0.2989 * input1.v(0)[i] + 0.5870 * input1.v(1)[i] + 0.1140 * input1.v(2)[i];
            image2[i] = 0.2989 * input2.v(0)[i] + 0.5870 * input2.v(1)[i] + 0.1140 * input2.v(2)[i];
        }
        
    } else
    {
        printf("ERROR :: src_ssim:: only grayscale or RGB images handled.\n");
        return EXIT_FAILURE;
    }
    
    
    // Read parameters
    float k1 = atof(oK1.value);
    float k2 = atof(oK2.value);
    int wsize = atoi(oSize.value);
    float wstd = atof(oStd.value);
    float L = atof(oRang.value);
    
    if((k1 < 0) || (k2 < 0))
    {
        printf("ERROR :: src_ssim:: k1 and k2 must be positive constant values.\n");
        return EXIT_FAILURE;
    }
    
    if(wsize < 2)
    {
        printf("ERROR :: src_ssim:: Gaussian kernel size must be grater than 1.\n");
        return EXIT_FAILURE;
        
    }
    
    if((width < wsize) || (height < wsize))
    {
        printf("ERROR :: src_ssim:: image width and height cannot be smaller than Gaussian kernel size.\n");
        return EXIT_FAILURE;
    }
    
    
    // Automatic downsampling
    int scale = MAX(1, roundf((float) (MIN(width, height) / 256.0f)));
    
    if(scale > 1)
    {
        float *aux1 = new float[dim];
        float *aux2 = new float[dim];
        
        libUSTG::fpCopy(image1, aux1, dim);
        libUSTG::fpCopy(image2, aux2, dim);
        
        delete[] image1;
        delete[] image2;
        
        int s_width = (int) floor((float) width / (float) scale);
        int s_height = (int) floor((float) height / (float) scale);
        int s_dim = s_width * s_height;
        
        image1 = new float[s_dim];
        image2 = new float[s_dim];
        
        libUSTG::fiImageSampleAglomeration(aux1, image1, scale, width, height);
        libUSTG::fiImageSampleAglomeration(aux2, image2, scale, width, height);
        
        width = s_width;
        height = s_height;
        dim = s_dim;
        
        delete[] aux1;
        delete[] aux2;
    }
    
    // SSIM constants
    float C1 = (k1 * L) * (k1 * L);
    float C2 = (k2 * L) * (k2 * L);
    
    
    // Normalized Gaussian window
    float *window = new float[wsize * wsize];
    libUSTG::fiFloatDirectionalGaussKernel(wstd, wstd, 0.0, window, wsize, wsize);
    
    
    // Compute image statistics
    float *mu1 = new float[dim];
    float *mu2 = new float[dim];
    
    libUSTG::fiConvol(image1, mu1, width, height, window, wsize, wsize, BOUNDARY_CONDITION_SYMMETRIC);
    libUSTG::fiConvol(image2, mu2, width, height, window, wsize, wsize, BOUNDARY_CONDITION_SYMMETRIC);
    
    float *mu1_sq = new float[dim];
    float *mu2_sq = new float[dim];
    float *mu12 = new float[dim];
    
    for(int i = 0; i < dim; i++)
    {
        mu1_sq[i] = mu1[i] * mu1[i];
        mu2_sq[i] = mu2[i] * mu2[i];
        mu12[i] = mu1[i] * mu2[i];
    }
    
    float *image1_sq = new float[dim];
    float *image2_sq = new float[dim];
    float *image12 = new float[dim];
    
    for(int i = 0; i < dim; i++)
    {
        image1_sq[i] = image1[i] * image1[i];
        image2_sq[i] = image2[i] * image2[i];
        image12[i] = image1[i] * image2[i];
    }
    
    float *sigma1_sq = new float[dim];
    float *sigma2_sq = new float[dim];
    float *sigma12 = new float[dim];
    
    
    libUSTG::fiConvol(image1_sq, sigma1_sq, width, height, window, wsize, wsize, BOUNDARY_CONDITION_SYMMETRIC);
    libUSTG::fiConvol(image2_sq, sigma2_sq, width, height, window, wsize, wsize, BOUNDARY_CONDITION_SYMMETRIC);
    libUSTG::fiConvol(image12, sigma12, width, height, window, wsize, wsize, BOUNDARY_CONDITION_SYMMETRIC);
    
    for(int i = 0; i < dim; i++)
    {
        sigma1_sq[i] = sigma1_sq[i] - mu1_sq[i];
        sigma2_sq[i] = sigma2_sq[i] - mu2_sq[i];
        sigma12[i] = sigma12[i] - mu12[i];
    }
    
    
    // Compute SSIM index
    float *ssim_map = new float[dim];
    
    if(C1 > fTiny && C2 > fTiny)
    {
        for(int i = 0; i < dim; i++)
            ssim_map[i] = ((2.0f * mu12[i] + C1) * (2.0f * sigma12[i] + C2)) / ((mu1_sq[i] + mu2_sq[i] + C1) * (sigma1_sq[i] + sigma2_sq[i] + C2));
    
    } else
    {
        libUSTG::fpClear(ssim_map, 1.0f, dim);
        
        for(int i = 0; i < dim; i++)
        {
            float num1 = 2.0f * mu12[i] + C1;
            float num2 = 2.0f * sigma12[i] + C2;
            float den1 = mu1_sq[i] + mu2_sq[i] + C1;
            float den2 = sigma1_sq[i] + sigma2_sq[i] + C2;
            
            if(den1 * den2 > fTiny)
                ssim_map[i] = (num1 * num2) / (den1 * den2);
            
            if(fabs(den1) > fTiny && fabs(den2) < fTiny)
                ssim_map[i] = num1 / den1;
        }
    }
    
    float ssim = libUSTG::fpMean(ssim_map, dim);
    
    printf("%1.4f\n", ssim);
    
    
    // Delete allocated memory
    delete[] image1;
    delete[] image2;
    delete[] window;
    delete[] mu1;
    delete[] mu2;
    delete[] mu1_sq;
    delete[] mu2_sq;
    delete[] mu12;
    delete[] image1_sq;
    delete[] image2_sq;
    delete[] image12;
    delete[] sigma1_sq;
    delete[] sigma2_sq;
    delete[] sigma12;
    delete[] ssim_map;
    
    
    return EXIT_SUCCESS;
}
