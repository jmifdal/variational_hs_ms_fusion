#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;



int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
	OptStruct oO = {"o", 0,  NULL, NULL,"flag for orthogonal transformation"};  options.push_back(&oO);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct py = {"y", NULL, "output image"}; parameters.push_back(&py);
	ParStruct pu = {"u", NULL, "output image"}; parameters.push_back(&pu);
	ParStruct pv = {"v", NULL, "output image"}; parameters.push_back(&pv);
	
	if (!parsecmdline("src_rgb_2_yuv", "Transforms RGB to YUV coordinates", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! input
	libUSTG::cflimage input;
	
	input.load(pinput.value);
	
	input.Rgb2Yuv(oO.flag);

	
	libUSTG::flimage y = input.getChannel(0);
	y.save( py.value);
	
	libUSTG::flimage u = input.getChannel(1);
	u.save( pu.value);
	
	libUSTG::flimage v = input.getChannel(2);
	v.save( pv.value);
	
	
	return 1;
	
	
}
