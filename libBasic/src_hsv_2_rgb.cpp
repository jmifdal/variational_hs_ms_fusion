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
    ParStruct ph = {"H", NULL, "input image"}; parameters.push_back(&ph);
    ParStruct ps = {"S", NULL, "input image"}; parameters.push_back(&ps);
    ParStruct pv = {"V", NULL, "input image"}; parameters.push_back(&pv);
    ParStruct pout = {"image", NULL, "output image"}; parameters.push_back(&pout);
    
    if(!parsecmdline("src_hsv_2_rgb", "Transforms HSV to RGB coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    //! Input
    libUSTG::flimage H, S, V;
    H.load(ph.value);
    S.load(ps.value);
    V.load(pv.value);
    
    
    //! Process
    libUSTG::cflimage output(H.w(), H.h(), 3);
    libUSTG::fiHsv2Rgb(H.v(), S.v(), V.v(), output.v(0), output.v(1), output.v(2), H.w() * H.h());
    
    
    //! Saving output
    output.save(pout.value);
    
    
    return EXIT_SUCCESS;
}