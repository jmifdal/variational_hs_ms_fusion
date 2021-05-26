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
    ParStruct pl = {"L", NULL, "output image"}; parameters.push_back(&pl);
    
    if(!parsecmdline("src_rgb_2_hsl", "Transforms RGB to HSL coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    //! Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    
    //! Process
    libUSTG::flimage H(input.w(), input.h());
    libUSTG::flimage S(input.w(), input.h());
    libUSTG::flimage L(input.w(), input.h());
    
    libUSTG::fiRgb2Hsl(input.v(0), input.v(1), input.v(2), H.v(), S.v(), L.v(), input.w() * input.h());
    
    
    //! Saving output
    H.save(ph.value);
    S.save(ps.value);
    L.save(pl.value);
    
    return EXIT_SUCCESS;
    
    
}