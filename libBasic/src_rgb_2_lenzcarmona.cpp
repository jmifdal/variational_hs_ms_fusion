#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
    
    std::vector<ParStruct *> parameters;
    ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
    ParStruct pl = {"l", NULL, "output L channel"}; parameters.push_back(&pl);
    ParStruct pc1 = {"c1", NULL, "output C1 channel"}; parameters.push_back(&pc1);
    ParStruct pc2 = {"c2", NULL, "output C2 channel"}; parameters.push_back(&pc2);
    
    if (!parsecmdline("src_rgb_2_lenzcarmona", "Transform from RGB to (L, C1, C2) - LenzCarmona coordinates", argc, argv,
                      options, parameters))
        return EXIT_FAILURE;
    
    
    //! Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    
    //! Process
    libUSTG::flimage L(input.w(), input.h());
    libUSTG::flimage C1(input.w(), input.h());
    libUSTG::flimage C2(input.w(), input.h());;
    
    libUSTG::fiRgb2LenzCarmona(input.v(0), input.v(1), input.v(2), L.v(), C1.v(), C2.v(), input.w() * input.h());
    
    
    //! Saving outputs
    L.save(pl.value);
    C1.save(pc1.value);
    C2.save(pc2.value);
    
    
    return EXIT_SUCCESS;
}