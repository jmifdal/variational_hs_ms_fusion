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
	ParStruct pinput = {"image", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	ParStruct px = {"x", NULL, "top left x coordinate"}; parameters.push_back(&px);
	ParStruct py = {"y", NULL, "top left y coordinate"}; parameters.push_back(&py);
	ParStruct pw = {"w", NULL, "width of extracted zone"}; parameters.push_back(&pw);
	ParStruct ph = {"h", NULL, "height of extracted zone"}; parameters.push_back(&ph);
	 
	
	if (!parsecmdline("src_extract"," extracts rectangular zone of an image", argc, argv, options, parameters))
		return 0;
	

	// parameters
	int x = atoi(px.value);
	int y = atoi(py.value);
	int w = atoi(pw.value);
	int h = atoi(ph.value);	

	// input
	libUSTG::cflimage input;
	input.load(pinput.value);
	

	// process
	libUSTG::cflimage output = input.copy(x, y, w, h);

	
	// save
	output.save(pout.value); 
	return 1;
	
}
