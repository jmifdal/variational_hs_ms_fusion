#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;



int main(int argc, char **argv)
{
	
    vector <OptStruct *> options;
    OptStruct oP = {"p", 0,  NULL, NULL, "modify one pixel"};  options.push_back(&oP);
	
	
	vector<ParStruct *> parameters;
	ParStruct pr = {"r", NULL, "input image"}; parameters.push_back(&pr);
	ParStruct pg = {"g", NULL, "input image"}; parameters.push_back(&pg);
	ParStruct pb = {"b", NULL, "input image"}; parameters.push_back(&pb);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	
	if (!parsecmdline("src_put_channels", "Puts rgb channels into color image", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! input
	libUSTG::flimage r,g,b;
	
	r.load(pr.value);
	g.load(pg.value);
	b.load(pb.value);
    
    int flagModify = oP.flag;
    
	
	//! output
	libUSTG::cflimage output(r,g,b);
    
    if(flagModify)
    {
        output.v(0)[0]=1.0f;
        output.v(1)[0]=255.0f;
        output.v(2)[0]=125.0f;
    }
    
	output.save(pout.value);
    	
	
	return EXIT_SUCCESS;
	
	
}
