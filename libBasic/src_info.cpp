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
	
	
	if (!parsecmdline("src_info", "prints image info", argc, argv, options, parameters))
		return 0;
	


	
//////////////////////////////////////////////// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	printf("input: %s \n", pinput.value);
	printf("width: %d\n", input.w());
	printf("height: %d\n", input.h());
	printf("channels: %d\n", input.c());
	
	
	for (int ii = 0; ii < input.c(); ii++)
	{
	
		printf("%d: (m,M) = (%f,%f)\n", ii, input.min_channel(ii), input.max_channel(ii) );
	}

	
	return 1;
}
