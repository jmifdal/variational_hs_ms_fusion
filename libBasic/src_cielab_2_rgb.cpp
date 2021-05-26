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
    ParStruct pl = {"L", NULL, "input image"}; parameters.push_back(&pl);
    ParStruct pa = {"a", NULL, "input image"}; parameters.push_back(&pa);
    ParStruct pb = {"b", NULL, "input image"}; parameters.push_back(&pb);
    ParStruct pout = {"image", NULL, "output image"}; parameters.push_back(&pout);
    
    if(!parsecmdline("src_cielab_2_rgb", "Transforms CIELab to RGB coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    //! Input
    libUSTG::flimage L, a, b;
    L.load(pl.value);
    a.load(pa.value);
    b.load(pb.value);
    
    //! Process
    libUSTG::cflimage output(L.w(), L.h(), 3);
    
    libUSTG::fiLab2Rgb(L.v(), a.v(), b.v(), output.v(0), output.v(1), output.v(2), L.w() * L.h());
  
    for(int c = 0; c < output.c(); c++)
        for(int i = 0; i < output.w() * output.h(); i++)
            output.v(c)[i] *= 255.0f;
    
    //! Saving output
    output.save(pout.value);
    
    
    return EXIT_SUCCESS;
}