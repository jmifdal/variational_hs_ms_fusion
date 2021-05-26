#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "../library/libImage.h"
#include "../library/libBasic.h"

using namespace std;


int main(int argc, char **argv)
{
	
	vector <OptStruct *> options;
    OptStruct oI = {"i:", 0,  NULL, NULL,"get only channel i"};  options.push_back(&oI);	
	
    
    
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct pout = {"output", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("src_get_channels", "gets channels of an image", argc, argv, options, parameters))
		return 0;
    
	
    
	// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
    if (!oI.flag)
    {
        //printf("... saving channels: \t"); fflush(stdout);
        
        for (int ii=0; ii < input.c(); ii++)
        {
            
            // image
            libUSTG::flimage finput(input.w(), input.h(), input.v(ii));
            
            
            // name of image
            stringstream s;
            s << ii;
            
            string nom(pout.value);
            nom += s.str();
            
            
            // save
            finput.save(nom.c_str());
            //printf("%s\t", nom.c_str());
            
        }
        
        //printf("\n");
        
	} else
    {
        
        int ii = atoi(oI.value);
        assert(ii >= 0 && ii < input.c());
        
        libUSTG::flimage finput(input.w(), input.h(), input.v(ii));
        finput.save(pout.value);
        //printf("... saving channel %d as image %s \n", ii, pout.value);
    }
	
	return 1;
	
}
