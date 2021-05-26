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
    ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
    ParStruct ph = {"H", NULL, "output image"}; parameters.push_back(&ph);
    ParStruct ps = {"S", NULL, "output image"}; parameters.push_back(&ps);
    ParStruct pv = {"V", NULL, "output image"}; parameters.push_back(&pv);
    
    if(!parsecmdline("src_rgb_2_hsv", "Transforms RGB to HSV coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    //! Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    
    //! Process
    libUSTG::flimage H(input.w(), input.h());
    libUSTG::flimage S(input.w(), input.h());
    libUSTG::flimage V(input.w(), input.h());
    
    libUSTG::fiRgb2Hsv(input.v(0), input.v(1), input.v(2), H.v(), S.v(), V.v(), input.w() * input.h());
    
    
    //! Saving output
    H.save(ph.value);
    S.save(ps.value);
    V.save(pv.value);
    
    return EXIT_SUCCESS;
    
    
}