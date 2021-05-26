#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#include "../library/libImage.h"
#include "../library/libBasic.h"


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


using namespace std;

int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
	OptStruct oL = {"l:", 0,  NULL, NULL,"rotate left"};  options.push_back(&oL);
	OptStruct oR = {"r:", 0,  NULL, NULL,"rotate right"};  options.push_back(&oR);
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	

	if (!parsecmdline("src_simple_transformation", "Rotate image left or right", argc, argv, options, parameters))
		return 0;


	
	/// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
    if (oL.flag)
    {
        libUSTG::cflimage output(input.h(), input.w(), input.c());
        for (int cc=0; cc < output.c(); cc++)
            for (int jj=0; jj < output.h(); jj++)
                for (int ii=0; ii < output.w(); ii++)
        {
            output[cc * output.c() + jj * output.w() + ii] = input[ cc * input.c() + ii * input.w() + input.w() - 1 - jj];
        }
        
        output.save(oL.value);
        
    }
	
    
    if (oR.flag)
    {
        libUSTG::cflimage output(input.h(), input.w(), input.c());
        for (int cc=0; cc < input.c(); cc++)
            for (int jj=0; jj < output.h(); jj++)
                for (int ii=0; ii < output.w(); ii++)
                {
                    output[cc * output.c() + jj * output.w() + ii] = input[ cc * input.c() + (input.h()-1-ii)* input.w() + jj];
                }
        
        output.save(oR.value);
        
    }
	
    
    
	return 1;
		
}
