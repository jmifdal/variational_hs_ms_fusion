#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;

int main(int argc, char* argv[])
{
    
    vector <OptStruct *> options;
    OptStruct oB = {"b:", 0, "0",   NULL, "boundary value if not interpolated"}; options.push_back(&oB);
    OptStruct oI  = {"i", 0, NULL,   NULL, "flag for using inverse homography"}; options.push_back(&oI);
    OptStruct oO  = {"o:", 0, "3",   NULL, "interpolation order"}; options.push_back(&oO);
    
    
    vector<ParStruct *> parameters;
    ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
    ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
    ParStruct pH = {"H", NULL, "homography"}; parameters.push_back(&pH);
    
    
    if (!parsecmdline("src_transform_homography" ,"Applies homography", argc, argv, options, parameters))
    return 0;
    
    
    
    //! Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    //! Homography
    std::ifstream ifilestr(pH.value);
    libUSTG::laMatrix H(3,3);
    
    if (!ifilestr.is_open())
    {
        printf("file %s impossible to open\n", pH.value);
        exit(-1);
    }
    
    
    for (int ii=0; ii < 3; ii++)
    for (int jj=0; jj < 3; jj++)
    {
        ifilestr >> H[ii][jj];
        
    }
    
    ifilestr.close();
    
    
    //! Inverse
    if (oI.flag)
    {
        
        libUSTG::laMatrix B(3,3);
        libUSTG::luinv(H,B);
        H=B;
    }
    
    
    //! Process
    libUSTG::cflimage output = input;
    
    float bg = atof(oB.value);
    int order = atoi(oO.value);
    
    //! normalize homography
    if (H[2][2] != 0.0f)
    for(int i=0; i < 3; i++)
    for(int j=0; j < 3; j++) H[i][j] /=  H[2][2];
    
    
    //! apply transormation tmp = H-1 (image)
    libUSTG::laMatrix Hinv = H;
    luinv(H, Hinv);
    
    
    //! Transform images
    if (order==3)
    {
        for (int ii=0; ii < input.c(); ii++)
        {
            bicubic_homography_interpolation(input.v(ii),input.w(),input.h(), Hinv, bg,  output.v(ii), output.w(), output.h());
        }
        
    } else
    {
        for (int ii=0; ii < input.c(); ii++)
        {
            apply_planar_homography(input.v(ii),input.w(),input.h(), Hinv, bg,order, output.v(ii), 0, 0, output.w(), output.h());
        }
    }
    
    
    //! Save
    output.save( pout.value);
    
}


