#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#include "../library/libImage.h"
#include "../library/libBasic.h"


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


using namespace std;

int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
	OptStruct oS = {"s:", 0,  NULL, NULL,"step value"};  options.push_back(&oS);	
	OptStruct oN = {"n:", 0,  NULL, NULL,"number of bins"};  options.push_back(&oN);	
	OptStruct oM = {"m:", 0,  NULL, NULL,"mask"};  options.push_back(&oM);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	

	if (!parsecmdline("src_histogram", "Image color histogram", argc, argv, options, parameters))
		return 0;


	
	/// Input
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	

	char flag;
	if (oS.flag && oN.flag) { printf("only one option -s or -n can be activated\n"); exit(-1);}

	if (!oS.flag && !oN.flag) { printf("please select one of -s or -n options\n"); exit(-1);}

	int iN;
	float fS;
	
	if (oS.flag) { fS = atof(oS.value); flag='s';}
	
	if (oN.flag) { iN = atoi(oN.value); flag='n';}
	
	
	float *mValues = NULL;
	int imaskSize = 0;
	libUSTG::flimage mask;
	
	if (oM.flag)
	{
		mask.load(oM.value);
		for (int ii=0; ii < mask.wh(); ii++) if (mask[ii] > 0.0f) imaskSize++;
		
		mValues = new float[imaskSize];
	}
	printf("imaskSize: %d\n", input.wh());
	printf("imaskSize: %d\n", imaskSize);
	
	
	for (int ii=0; ii < input.c(); ii++)
	{
	
		
		
        libUSTG::laVector histo = NULL;
        
		if (oM.flag)
		{
			
			float *iptr = input.v(ii);
			
			int kk=0;
			for (int jj=0; jj < input.wh(); jj++)
			if (mask[jj] > 0.0)
			{
				mValues[kk] = iptr[jj];
				kk++;
			}
			
			
			libUSTG::fpHisto(mValues, histo, NULL, NULL, &iN, &fS, imaskSize, flag);
		
		}
		else
			libUSTG::fpHisto(input.v(ii), histo, NULL, NULL, &iN, &fS, input.wh(), flag);
		
		
		
		
		// name of histogram file
		stringstream s;
		s << ii;
		
		string nom(pout.value);
		nom += s.str();
		
		
		// save
		ofstream ofilestr(nom.c_str());
		
		if (!ofilestr.is_open())
		{
			printf("file %s impossible to open\n", nom.c_str());
			exit(-1); 
		}
		
		
		float fX = input.min_channel(ii) + 0.5f * fS ;
		for(int ij = 0; ij < iN; ij++)
		{
			ofilestr << fX  << " " << histo[ij] << endl;
			fX += fS;
		}
		
		
		histo.erase();
		
	}
	
	
	
	return 1;
		
}
