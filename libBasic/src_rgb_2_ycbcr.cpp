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
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct py = {"y", NULL, "output image"}; parameters.push_back(&py);
    ParStruct pcb = {"cb", NULL, "output image"}; parameters.push_back(&pcb);
    ParStruct pcr = {"cr", NULL, "output image"}; parameters.push_back(&pcr);
	
	if (!parsecmdline("fip_rgb_2_yuv","Transform to yuv coordinates", argc, argv, options, parameters))
		return 0;
	
	
	
	
	//! input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	libUSTG::flimage y(input.w(), input.h());
	libUSTG::flimage cb(input.w(), input.h());
	libUSTG::flimage cr(input.w(), input.h());

    
    /* Conversion for values in (0,1)
     origT = [65.481 128.553 24.966;...
             -37.797 -74.203 112; ...
             112 -93.786 -18.214];
    origOffset = [16;128;128];
    ycbcr = origT * rgb + origOffset
     */
    
    for (int i=0; i < input.wh(); i++)
    {
        
        float r = input[i];
        float g = input[input.wh() + i];
        float b = input[2 * input.wh() + i];
        
        r /= 255.0f; if (r<0.0f) r=0.0f; if (r>1.0f) r=1.0f;
        g /= 255.0f; if (g<0.0f) g=0.0f; if (g>1.0f) g=1.0f;
        b /= 255.0f; if (b<0.0f) b=0.0f; if (b>1.0f) b=1.0f;
        
        y[i] = 65.481 * r + 128.553 * g + 24.966 * b + 16.0f;
        cb[i] = -37.797 * r - 74.203 * g + 112 * b + 128.0f;
        cr[i] = 112 * r - 93.786 * g - 18.214 * b + 128.0f;
        
        
        
    }
    
    
    
    y.save( py.value);
    cb.save( pcb.value);
    cr.save( pcr.value);
	
	
	return 1;
	
	
}
