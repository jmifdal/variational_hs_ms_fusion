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
    ParStruct pinput = {"image", NULL, "input image"}; parameters.push_back(&pinput);
    ParStruct pl = {"L", NULL, "output image"}; parameters.push_back(&pl);
    ParStruct pa = {"a", NULL, "output image"}; parameters.push_back(&pa);
    ParStruct pb = {"b", NULL, "output image"}; parameters.push_back(&pb);
    
    if(!parsecmdline("src_rgb_2_cielab", "Transforms RGB to CIELab coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    //! Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    for(int c = 0; c < input.c(); c++)
        for(int i = 0; i < input.w() * input.h(); i++)
            input.v(c)[i] /= 255.0f;
    
    
    //! Process
    libUSTG::flimage L(input.w(), input.h());
    libUSTG::flimage a(input.w(), input.h());
    libUSTG::flimage b(input.w(), input.h());
    
    libUSTG::fiRgb2Lab(input.v(0), input.v(1), input.v(2), L.v(), a.v(), b.v(), input.w() * input.h());
    
    
    //! Saving outputs
    L.save(pl.value);
    a.save(pa.value);
    b.save(pb.value);

    
    
    
    return EXIT_SUCCESS;
}