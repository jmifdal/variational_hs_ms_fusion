#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
	OptStruct oX = {"x:", 0,  NULL, NULL,"x gradient"};  options.push_back(&oX);	
	OptStruct oY = {"y:", 0,  NULL, NULL,"y gradient"};  options.push_back(&oY);	

    
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pgrad = {"gradient", NULL, "output file"}; parameters.push_back(&pgrad);
	ParStruct pori = {"orientation", NULL, "threshold"}; parameters.push_back(&pori);
	
	
	if (!parsecmdline("src_gradient", "computes image gradient and orientation", argc, argv, options, parameters))
		return 0;
	
	
	
	
	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	//! Process
	libUSTG::cflimage orientation(input.w(), input.h(), input.c());
	libUSTG::cflimage gradx(input.w(), input.h(), input.c());
	libUSTG::cflimage grady(input.w(), input.h(), input.c());
    
    
    libUSTG::cflimage gradient = input.gradient(gradx, grady, orientation, 'f');
	
	
	//! Save
	gradient.save(pgrad.value);
	orientation.save(pori.value);

    if (oX.flag) gradx.save(oX.value);
    if (oY.flag) grady.save(oY.value);
    


}

