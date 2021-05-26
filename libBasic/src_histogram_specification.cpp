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
    ParStruct pinput1 = {"reference", NULL, "reference input image"}; parameters.push_back(&pinput1);
    ParStruct pinput2 = {"modified", NULL, "input image to be modified"}; parameters.push_back(&pinput2);
    ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
    
    if(!parsecmdline("src_histogram_specification", "Global histogram specification",
                      argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // input
    libUSTG::cflimage input1;
    input1.load(pinput1.value);
    
    
    // input
    libUSTG::cflimage input2;
    input2.load(pinput2.value);
    
    
    // tests
    assert(input1.c() == input2.c());
    
    
    // process
    for (int i=0; i < input1.c(); i++)
        libUSTG::fk_histogram_specification(input1.v(i), input2.v(i), input2.v(i), input1.w(), input1.h(), input2.w(), input2.h());
    
    input2.save(pout.value);
    

    return EXIT_SUCCESS;
}
