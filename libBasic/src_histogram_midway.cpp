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
	ParStruct pinput1 = {"image", NULL, "input image1"}; parameters.push_back(&pinput1);
	ParStruct pinput2 = {"image", NULL, "input image2"}; parameters.push_back(&pinput2);
	ParStruct pout1 = {"out", NULL, "output image1"}; parameters.push_back(&pout1);
	ParStruct pout2 = {"out", NULL, "output image2"}; parameters.push_back(&pout2);
	
	
	if (!parsecmdline("src_histogram_midway", "histogram midway equalization", argc, argv, options, parameters))
		return 0;
	
	
	
	
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
	{
	
        libUSTG::fk_histogram_midway(input1.v(i), input2.v(i), input1.v(i), input2.v(i), input1.w(), input1.h(), input2.w(), input2.h());
		
	}
	
	
	input1.save(pout1.value);
	input2.save(pout2.value);
	



}
