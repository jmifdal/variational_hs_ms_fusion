

#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;


int main(int argc, char **argv)
{

	vector<OptStruct *> options;
	OptStruct ocr = {"r:", 0, "0.333333", NULL,"r weight value"};  options.push_back(&ocr);	
	OptStruct ocg = {"g:", 0, "0.333333", NULL,"g weight value"};  options.push_back(&ocg);	
	OptStruct ocb = {"b:", 0, "0.333333", NULL,"b weight value"};  options.push_back(&ocb);
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image1", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);

	
	if (!parsecmdline("src_gray", "Converts color image to gray", argc, argv, options, parameters))
		return 0;


	
	//! Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	//!  Process
	libUSTG::flimage output = input.getGray(atof(ocr.value), atof(ocg.value), atof(ocb.value));


	//! Save
	output.save( pout.value); 
	return 1;
	
}
