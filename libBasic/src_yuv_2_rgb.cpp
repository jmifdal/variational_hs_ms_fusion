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
	ParStruct py = {"y", NULL, "input image"}; parameters.push_back(&py);
	ParStruct pu = {"u", NULL, "input image"}; parameters.push_back(&pu);
	ParStruct pv = {"v", NULL, "input image"}; parameters.push_back(&pv);
	ParStruct poutput = {"image", NULL, "output image"}; parameters.push_back(&poutput);
	
	if (!parsecmdline("src_yuv_2_rgb", "Transforms from YUV to RGB coordinates", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! input
	libUSTG::flimage y;
	y.load(py.value);

	libUSTG::flimage u;
	u.load(pu.value);

	libUSTG::flimage v;
	v.load(pv.value);
	
	
	//! process
	libUSTG::cflimage yuv(y,u,v);
	yuv.Yuv2Rgb(oO.flag);
	

	
	yuv.save(poutput.value);
	
	return 1;
	
}
