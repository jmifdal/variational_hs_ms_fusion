#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char **argv)
{
    
	vector <OptStruct *> options;
	OptStruct oP = {"p:", 0,  "2", NULL,"L^p distance (p=0,1,2)"};  options.push_back(&oP);
	OptStruct oM = {"m:", 0,  NULL, NULL, "mask>0 points used in distance"};  options.push_back(&oM);
	OptStruct oV = {"v", 0,  NULL, NULL, "verbose mode"};  options.push_back(&oV);
	OptStruct oB = {"b:", 0,  "0", NULL, "boundary elimination "};  options.push_back(&oB);
	OptStruct oC = {"c:", 0,  "0", NULL, "distance for specific channels"};  options.push_back(&oC);
	
    
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image", NULL, "image"}; parameters.push_back(&pinput2);
	
	if (!parsecmdline("src_dist_lp", "computes distance l^p", argc, argv, options, parameters))
		return 0;
	
	
    
	// Parameters
	int p = atoi(oP.value);
    
    
	// Input
    libUSTG::cflimage input;
	input.load(pinput.value);
	
    
	
	// Input
	libUSTG::cflimage input2;
	input2.load(pinput2.value);
	
    
	
	// Mask
	libUSTG::flimage imask;
	if (oM.flag)
		imask.load(oM.value);
	
    
    
	// process
	assert(input.c() == input2.c() && input.w() == input2.w() && input.h() == input2.h());
	
	
	// remove boundary
	int boundary = atoi(oB.value);
    
	libUSTG::cflimage cinput = input.copy(boundary/2, boundary/2, input.w() - boundary, input.h()-boundary);
	libUSTG::cflimage cinput2 = input2.copy(boundary/2, boundary/2, input2.w() - boundary, input2.h()-boundary);
	
    
    
    
	libUSTG::flimage cimask;
	if (oM.flag)
		cimask = imask.copy(boundary/2, boundary/2, imask.w() - boundary, imask.h()-boundary);
    
	
    if (!oC.flag)
    {
        
        float fDist = 0.0f;
        for (int ii=0; ii < input.c(); ii++)
        {
            
            float fDif = 0.0f;
            
            if (oM.flag)
                fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), cimask.v(), p, cinput.wh());
            
            else
                fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), p, cinput.wh());
            
            fDist += fDif;
            
            
            if (oV.flag)
                printf("channel %d: %2.5f\n", ii, fDif);
        }
        
        
        fDist /= (float) cinput.c();
        
        
        if (oV.flag)
            printf("mean error: %2.5f\n",fDist);
        else
            printf("%2.5f\n", fDist);
    } else
    {
        
        int ii = atoi(oC.value);
        assert(ii>=0 && ii<input.c());
        
        float fDif = 0.0f;
        if (oM.flag)
			fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), cimask.v(), p, cinput.wh());
		
		else
			fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), p, cinput.wh());
		
        printf("channel %d: %2.5f\n", ii, fDif);
        
    }
    
	
	return 1;
	
	
}
