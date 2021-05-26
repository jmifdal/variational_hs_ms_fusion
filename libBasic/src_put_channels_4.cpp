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
	ParStruct pr = {"r", NULL, "input image"}; parameters.push_back(&pr);
	ParStruct pg = {"g", NULL, "input image"}; parameters.push_back(&pg);
    ParStruct pg2 = {"g2", NULL, "input image"}; parameters.push_back(&pg2);
    ParStruct pb = {"b", NULL, "input image"}; parameters.push_back(&pb);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	
	if (!parsecmdline("src_put_channels_4", "Put rgb-ninfrared channels into color image", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! input
	libUSTG::flimage r,g,g2,b;
	
	r.load(pr.value);
	g.load(pg.value);
    g2.load(pg2.value);
	b.load(pb.value);
	
	//! input
	libUSTG::cflimage output(r.w(), r.h(), 4);
    

    libUSTG::fpCopy(r.v(), output.v(0), output.wh());
    libUSTG::fpCopy(g.v(), output.v(1), output.wh());
    libUSTG::fpCopy(g2.v(), output.v(2), output.wh());
    libUSTG::fpCopy(b.v(), output.v(3), output.wh());
    
        
	output.save(pout.value);
    	
	
	return 1;
	
	
}
