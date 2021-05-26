#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libBasic.h"
#include "../library/libImage.h"


//BEWARE: we fix the first input image to be the reference one

int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
	OptStruct oB = {"b:", 0,  "0", NULL, "boundary elimination "};  options.push_back(&oB);
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image", NULL, "image"}; parameters.push_back(&pinput2);
	
	if(!parsecmdline("src_psnr","psnr", argc, argv, options, parameters))
		return EXIT_FAILURE;
	
	
	//! input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	
	//! input2
	libUSTG::cflimage input2;
	input2.load(pinput2.value);
    
    assert(input.c() == input2.c() && input.w() == input2.w() && input.h() == input2.h());

	
    // remove boundary
	int boundary = atoi(oB.value);
	libUSTG::cflimage cinput = input.copy(boundary/2, boundary/2, input.w() - boundary, input.h()-boundary);
	libUSTG::cflimage cinput2 = input2.copy(boundary/2, boundary/2, input2.w() - boundary, input2.h()-boundary);
	
    //! process
    float fDist = 0.0f;
    float fPSNR = 0.0f;

    //Recovering the maximum
    float max_ref=cinput.max();
    float max_ref_2=max_ref*max_ref;
    float ratio=0.0f;

	
	for (int jj=0; jj < cinput.whc(); jj++)
	{
		
		float dif = cinput[ jj] - cinput2[ jj];
		fDist += dif * dif;
		
	}
	
    fDist /= (float) (cinput.whc());
    //fDist /= sqrtf(fDist);
	
    ratio=(float)max_ref_2/fDist;

    //fPSNR = 10.0f * log10f(255.0f * 255.0f / (fDist * fDist));
    fPSNR = 10.0f * log10f(ratio);
	
    //printf("PSNR: %2.2f\n", fPSNR);
    printf("%2.5f\n", fPSNR);
	
    
    return EXIT_SUCCESS;
}





/****************************************************
 * sauvegarde du code original : 11/06/2019 (Jamila)
 * **************************************************/

/*#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
	OptStruct oB = {"b:", 0,  "0", NULL, "boundary elimination "};  options.push_back(&oB);

    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image", NULL, "image"}; parameters.push_back(&pinput2);

	if(!parsecmdline("src_psnr","psnr", argc, argv, options, parameters))
		return EXIT_FAILURE;


	//! input
	libUSTG::cflimage input;
	input.load(pinput.value);


	//! input2
	libUSTG::cflimage input2;
	input2.load(pinput2.value);

    assert(input.c() == input2.c() && input.w() == input2.w() && input.h() == input2.h());


    // remove boundary
	int boundary = atoi(oB.value);
	libUSTG::cflimage cinput = input.copy(boundary/2, boundary/2, input.w() - boundary, input.h()-boundary);
	libUSTG::cflimage cinput2 = input2.copy(boundary/2, boundary/2, input2.w() - boundary, input2.h()-boundary);

    //! process
    float fDist = 0.0f;
    float fPSNR = 0.0f;

	for (int jj=0; jj < cinput.whc(); jj++)
	{

		float dif = cinput[ jj] - cinput2[ jj];
		fDist += dif * dif;

	}



    fDist /= (float) (cinput.whc());
    fDist /= sqrtf(fDist);

    fPSNR = 10.0f * log10f(255.0f * 255.0f / (fDist * fDist));

//    printf("PSNR: %2.2f\n", fPSNR);
    printf("%2.2f\n", fPSNR);


    return EXIT_SUCCESS;
}

 */
