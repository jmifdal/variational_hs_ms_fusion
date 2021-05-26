#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
    
    std::vector<ParStruct *> parameters;
    ParStruct pl = {"l", NULL, "input Y channel"}; parameters.push_back(&pl);
    ParStruct pc1 = {"c1", NULL, "input U channel"}; parameters.push_back(&pc1);
    ParStruct pc2 = {"c2", NULL, "input V channel"}; parameters.push_back(&pc2);
    ParStruct poutput = {"output", NULL, "output RGB image"}; parameters.push_back(&poutput);
    
    if(!parsecmdline("src_conrat_2_rgb", "Transform from (L,C1,C2) - Conrat to RGB coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    //! Inputs
    libUSTG::flimage L, C1, C2;
    L.load(pl.value);
    C1.load(pc1.value);
    C2.load(pc2.value);
    
    
    //! Process
    libUSTG::cflimage output(L.w(), L.h(), 3);
    libUSTG::fiConrat2Rgb(output.v(0), output.v(1), output.v(2), L.v(), C1.v(), C2.v(), L.w()*L.h());
    
    
    //! Saving output
    output.save(poutput.value);
    
    
    return EXIT_SUCCESS;
}