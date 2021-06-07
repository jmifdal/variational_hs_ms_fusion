#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "../library/libImage.h"
#include "../library/libHyperspectralWeights.h"
//#include "./libHyperspectralWeights.h"

int main(int argc, char **argv)
{
	std::vector <OptStruct *> options;
	OptStruct oHPatchDist = {"c:", 0, "2.5", NULL, "filtering parameter for patch distance in NL weights"}; options.push_back(&oHPatchDist);
	OptStruct oHSpatDist = {"d:", 0, "2.5", NULL, "filtering parameter for spatial distance in NL weights"}; options.push_back(&oHSpatDist);
	OptStruct oNneig = {"n:", 0, "25", NULL, "number of (most similar) neighboring pixels used in NL weights"}; options.push_back(&oNneig);
	OptStruct oBandWindSize = {"b:", 0, "1", NULL, "half-size of spectral window in NL weights"}; options.push_back(&oBandWindSize);
	OptStruct oResWindSize = {"r:", 0, "7", NULL, "half-size of research window in NL weights"}; options.push_back(&oResWindSize);
	OptStruct oPatchSize = {"p:", 0, "1", NULL, "half-size of patches in NL weights"}; options.push_back(&oPatchSize);

	std::vector<ParStruct *> parameters;
	ParStruct pMulti = {"multispectral", NULL, "input high-resolution MS image"};
	parameters.push_back(&pMulti);
	ParStruct pHyperInt = {"hyperspectral_int", NULL, "input spatially upsampled HS image"};
	parameters.push_back(&pHyperInt);
	ParStruct pSmatrix = {"Smatrix", NULL, "input spectral-downsampling matrix operator"};
	parameters.push_back(&pSmatrix);
	ParStruct pWxy = {"wxy", NULL, "output wxy NL weights"};
	parameters.push_back(&pWxy);
	ParStruct pWyx = {"wyx", NULL, "output wxy NL weights"};
	parameters.push_back(&pWyx);
	ParStruct pPosxy = {"posxy", NULL, "output posxy NL weights"};
	parameters.push_back(&pPosxy);
	ParStruct pPosyx = {"posyx", NULL, "output posyx NL weights"};
	parameters.push_back(&pPosyx);
	ParStruct pPosw = {"posw", NULL, "output posw NL weights"};
	parameters.push_back(&pPosw);
	ParStruct pNumNeighxy = {"numNeighbours", NULL, "output number of neighbours per pixel in wxy weights"};
	parameters.push_back(&pNumNeighxy);
	ParStruct pNumNeighyx = {"numNeighbours", NULL, "output number of neighbours per pixel in wyx weights"};
	parameters.push_back(&pNumNeighyx);

    
	char* arg = "INFO: Computation of the nonlocal weights.\n";

	if(!parsecmdline("src_hyperspectral_fusion_L2_weights", arg, argc, argv, options, parameters))
		return EXIT_FAILURE;


	// Parameters
	float hPatch = atof(oHPatchDist.value);
	float hSpatial = atof(oHSpatDist.value);
	int numNeighbours = atoi(oNneig.value);
	int bandwind = atoi(oBandWindSize.value);
	int reswind = atoi(oResWindSize.value);
	int compwind = atoi(oPatchSize.value);

	// Image data
	libUSTG::cflimage mimg, himgInt;
	mimg.load(pMulti.value);
	himgInt.load(pHyperInt.value);

	int width = mimg.w();
	int height = mimg.h();
	int dim = width * height;
	int ms_channels = mimg.c();
	int hs_channels = himgInt.c();

	if((himgInt.w() != width) || (himgInt.h() != height))
	{
		printf("ERROR :: src_hyperspectral_fusion_L2_weights :: uspampled HS image must be defined on high-resolution domain.\n");
		return EXIT_FAILURE;
	}

	float **fM = new float*[ms_channels];
	for(int m = 0; m < ms_channels; m++)
	{
		fM[m] = new float[dim];
		libUSTG::fpCopy(mimg.v(m), fM[m], dim);
	}

	float **fHint = new float*[hs_channels];
	for(int h = 0; h < hs_channels; h++)
	{
		fHint[h] = new float[dim];
		libUSTG::fpCopy(himgInt.v(h), fHint[h], dim);
	}

	// Spectral-downsampling operator
	libUSTG::flimage Simg;
	Simg.load(pSmatrix.value);

	if(Simg.c() != 1)
	{
		printf("ERROR :: src_hyperspectral_fusion_L2_weights :: spectral-downsampling matrix must be a one-channel image.\n");
		return EXIT_FAILURE;
	}

	if(Simg.w() != hs_channels)
	{
		printf("ERROR :: src_hyperspectral_fusion_L2_weights :: spectral-downsampling matrix width must be the number of HS channels.\n");
		return EXIT_FAILURE;
	}

	if(Simg.h() != ms_channels)
	{
		printf("ERROR :: src_hyperspectral_fusion_L2_weights :: spectral-downsampling matrix height must be the number of MS channels.\n");
		return EXIT_FAILURE;
	}

	float **S = new float*[ms_channels];

	for(int m = 0; m < ms_channels; m++)
	{
		S[m] = new float[hs_channels];

		for(int h = 0; h < hs_channels; h++)
			S[m][h] = Simg[m * hs_channels + h];
	}


	// Compute and save NL weights and vector with number of neighbours per pixel
	int flagNormalize = 1;
	int *numNeighxy = new int[dim];
	int *numNeighyx = new int[dim];


    std::vector< std::vector <float> > wxy;
    std::vector< std::vector <int> > posxy;
    std::vector< std::vector <float> > wyx;
    std::vector< std::vector <int> > posyx;
    std::vector< std::vector <int> > posw;

    libUSTGHYPER:: nlweights_ms_Smatrix_new(mimg, S, wxy, posxy, wyx, posyx, posw, hs_channels, hPatch,
            hSpatial, numNeighbours, reswind, compwind, flagNormalize);


    libUSTG::cflimage outWxy(numNeighbours, dim, hs_channels);
    libUSTG::cflimage outWyx(numNeighbours, dim, hs_channels);
    libUSTG::cflimage outPosxy(numNeighbours, dim, hs_channels);
    libUSTG::cflimage outPosyx(numNeighbours, dim, hs_channels);
    libUSTG::cflimage outPosw(numNeighbours, dim, hs_channels);

    for(int r = 0; r < dim; r++)
    {
        numNeighxy[r] = (int) wxy[r].size();
        numNeighyx[r] = (int) wyx[r].size();
    }

    for(int r = 0; r < hs_channels * dim - 1; r++)
    {
        int h = (int) (r / dim);
        int i = (int) (r % dim);

        for(int w = 0; w < (int) wxy[r].size(); w++)
        {
            int index = i * numNeighbours + w;
            outWxy.v(h)[index] = wxy[r][w];
            outPosxy.v(h)[index] = posxy[r][w];
        }

        for(int w = 0; w < (int) wyx[r].size(); w++)
        {
            int index = i * numNeighbours + w;
            outWyx.v(h)[index] = wyx[r][w];
            outPosyx.v(h)[index] = posyx[r][w];
            outPosw.v(h)[index] = posw[r][w];
        }
    }

    outWxy.save(pWxy.value);
    outWyx.save(pWyx.value);
    outPosxy.save(pPosxy.value);
    outPosyx.save(pPosyx.value);
    outPosw.save(pPosw.value);



	// Save images with number of neighbours per pixel
	libUSTG::flimage imnumxy(width, height);
	libUSTG::flimage imnumyx(width, height);

	for(int i = 0; i < dim; i++)
	{
		imnumxy.v()[i] = numNeighxy[i];
		imnumyx.v()[i] = numNeighyx[i];
	}

	imnumxy.save(pNumNeighxy.value);
	imnumyx.save(pNumNeighyx.value);

	// Delete allocated memory
	for(int m = 0; m < ms_channels; m++)
	{
		delete[] fM[m];
		delete[] S[m];
	}

	for(int h = 0; h < hs_channels; h++)
		delete[] fHint[h];

	delete[] fM;
	delete[] fHint;
	delete[] S;
	delete[] numNeighxy;
	delete[] numNeighyx;

	return EXIT_SUCCESS;
}
