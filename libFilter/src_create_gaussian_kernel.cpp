#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"



using namespace std;

int main(int argc, char **argv)
{
	
	
    vector <OptStruct *> options;
	OptStruct ox = {"x:", 0, "1.0", NULL, "gaussian standard deviation in x direction"};  options.push_back(&ox);	
	OptStruct oy = {"y:", 0, "1.0", NULL, "gaussian standard deviation in y direction"};  options.push_back(&oy);	
	OptStruct oa = {"a:", 0, "0.0", NULL, "filter angle"};  options.push_back(&oa);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pout = {"out", NULL, "output filter"}; parameters.push_back(&pout);
	
	if (!parsecmdline("src_create_gaussian_kernel", "Creates directional gaussian filter",
                      argc, argv, options, parameters))
		return 0;
	
	

    //! Parameters
 	float xsigma = atof(ox.value);
 	float ysigma = atof(oy.value);
	float angle =  atof(oa.value);
	
	int ksize = 4 * (int) ceilf(MAX(xsigma, ysigma)) + 1;
	
	
	int kwidth = ksize, kheight = ksize;
	
	
    libUSTG::flimage kernel(kwidth, kheight);
	
	
	libUSTG::fiFloatDirectionalGaussKernel(xsigma,ysigma,angle, kernel.v(0), kwidth, kheight);


	kernel.save(pout.value);
	
}
