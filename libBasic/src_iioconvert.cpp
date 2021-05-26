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
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);

	
	if (!parsecmdline("src_iioconvert", "Converts image format", argc, argv, options, parameters))
		return 0;


	
	//! Read input
    int d_w, d_h, d_c;
    float *d_v;
    
    d_v = iio_read_image_float_split(pinput.value, &d_w, &d_h, &d_c);
    
    libUSTG::flimage output(d_w, d_h, d_v);
    output.save(pout.value);
    
    return 1;
	
}
