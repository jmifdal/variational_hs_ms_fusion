#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char **argv)
{
    vector <OptStruct *> options;
    OptStruct oM = {"m:", 0,  NULL, NULL, "mask>0 points used in distance"};  options.push_back(&oM);
    OptStruct oB = {"b:", 0,  "0", NULL, "boundary elimination "};  options.push_back(&oB);
    
    vector<ParStruct *> parameters;
    ParStruct pinput1 = {"image1", NULL, "image1"}; parameters.push_back(&pinput1);
    ParStruct pinput2 = {"image2", NULL, "image2"}; parameters.push_back(&pinput2);
    
    if(!parsecmdline("src_sam", "compute SAM index", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Inputs
    libUSTG::cflimage input1, input2;
    input1.load(pinput1.value);
    input2.load(pinput2.value);
    
    assert(input1.c() == input2.c() && input1.w() == input2.w() && input1.h() == input2.h());
    
    
    // Mask
    libUSTG::flimage imask;
    
    if(oM.flag)
        imask.load(oM.value);
    
    
    // Remove boundary
    int boundary = atoi(oB.value);
    
    libUSTG::cflimage cinput1 = input1.copy(boundary/2, boundary/2, input1.w() - boundary, input1.h()-boundary);
    libUSTG::cflimage cinput2 = input2.copy(boundary/2, boundary/2, input2.w() - boundary, input2.h()-boundary);
    
    
    libUSTG::flimage cimask;
    if(oM.flag)
        cimask = imask.copy(boundary/2, boundary/2, imask.w() - boundary, imask.h()-boundary);
    
    int width = cinput1.w();
    int height = cinput1.h();
    int dim = width * height;
    int num_channels = cinput1.c();
    
    float **im1 = new float*[num_channels];
    float **im2 = new float*[num_channels];
    float *mask = new float[dim];
    
    for(int c = 0; c < num_channels; c++)
    {
        im1[c] = new float[dim];
        libUSTG::fpCopy(cinput1.v(c), im1[c], dim);
        
        im2[c] = new float[dim];
        libUSTG::fpCopy(cinput2.v(c), im2[c], dim);
    }
    
    if(oM.flag)
        libUSTG::fpCopy(cimask.v(), mask, dim);
    
    
    // Compute SAM value
    float sam = 0.0f;
    int num_pixels = 0;
    
    for(int i = 0; i < dim; i++)
    {
        if((oM.flag && mask[i] > 0) || (!oM.flag))
        {
            float norm1 = 0.0f;
            float norm2 = 0.0f;
            float prodscalar = 0.0f;
            
            for(int c = 0; c < num_channels; c++)
            {
                float im1val = im1[c][i];
                float im2val = im2[c][i];
                
                norm1 += (im1val * im1val);
                norm2 += (im2val * im2val);
                prodscalar += (im1val * im2val);
            }
            
            norm1 = sqrtf(norm1);
            norm2 = sqrtf(norm2);
            
            float div = norm1 * norm2;
            
            if(div > fTiny)
            {
                float arg = MIN(prodscalar / div, 1.0f);
                sam += acosf(arg);
                num_pixels++;
            }
        }
    }
    
    if(num_pixels > 0)
        sam /= (float) num_pixels;
    
    sam *= (180.0f / PI);
    
    printf("%2.5f\n", sam);
    
    
    // Free memory
    for(int c = 0; c < num_channels; c++)
    {
        delete[] im1[c];
        delete[] im2[c];
    }
    
    delete[] im1;
    delete[] im2;
    delete[] mask;
    
    return EXIT_SUCCESS;
}
