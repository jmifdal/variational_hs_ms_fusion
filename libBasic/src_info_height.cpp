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
	ParStruct pinput = {"image1", NULL, "image"}; parameters.push_back(&pinput);
	
	
	if (!parsecmdline("src_info_height","print image info", argc, argv, options, parameters))
		return 0;
	


	
//////////////////////////////////////////////// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	printf("%d\n", input.h());
		
	return 1;
}
