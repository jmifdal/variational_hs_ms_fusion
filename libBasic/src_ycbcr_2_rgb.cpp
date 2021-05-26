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
	ParStruct py = {"y", NULL, "input image"}; parameters.push_back(&py);
    ParStruct pcb = {"cb", NULL, "input image"}; parameters.push_back(&pcb);
    ParStruct pcr = {"cr", NULL, "input image"}; parameters.push_back(&pcr);
    ParStruct pOutput = {"image", NULL, "image"}; parameters.push_back(&pOutput);
	
	if(!parsecmdline("src_ycbcr_2_rgb","Transform to yuv coordinates", argc, argv, options, parameters))
        return EXIT_FAILURE;
	
	
	
	
	//! input
    libUSTG::flimage yim, cbim, crim;
    yim.load(py.value);
    cbim.load(pcb.value);
    crim.load(pcr.value);
    
    int dim = yim.wh();
    
    float *r = new float[dim];
    float *g = new float[dim];
    float *b = new float[dim];
    
    for (int i=0; i < dim; i++)
    {
        
        float y = yim[i];
        float cb = cbim[i];
        float cr = crim[i];
        
        y = 255.0f * (y - 16.0f);
        cb = 255.0f * (cb - 128.0f);
        cr = 255.0f * (cr - 128.0f);
        
        
        
       // r[i] = 0.0046f * y +  1.28e-09 * cb + 0.0063 * cr;
       // g[i] = 0.0046f * y - 0.0015 * cb - 0.0032 * cr;
       // b[i] = 0.0046f * y + 0.0079 * cb - 1.1977e-08 * cr;
        
        r[i] = 0.0046f * y  + 0.0063 * cr;
        g[i] = 0.0046f * y - 0.0015 * cb - 0.0032 * cr;
        b[i] = 0.0046f * y + 0.0079 * cb;
    }
    
    
    libUSTG::cflimage output(yim.w(), yim.h(), r, g, b);
    output.save(pOutput.value);
    
	
    delete[] r; delete[] g; delete[] b;
    
    
	return EXIT_SUCCESS;
	
	
}
