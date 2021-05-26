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
	ParStruct oNLwOpt = {"weightOpt", NULL, "input NL weights option\n"};
	parameters.push_back(&oNLwOpt);
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

	char* arg = "INFO: Computation of the nonlocal weights.\n"
			"      Options available:\n"
			"      - 1 : weights computed on MS (the same for all bands).\n"
			"      - 2 : weights computed on MS with Smatrix coefficients.\n"
			"      - 3 : weights computed on each upsampled HS band.\n"
			"      - 4 : weights computed on full upsampled HS (the same for all bands).\n"
			"      - 5 : weights computed on nearest upsampled HS bands.\n"
			"      - 6 : neighbors selected on full upsampled HS and weights computed just on associated HS band.\n"
			"      - 7 : neighbors selected on full MS and weights computed just on associated HS band.\n"
			"      - 8 : neighbors selected on full MS with Smatrix coefficients and weights computed just on associated HS band.\n"
			"	   - 9 : weights computed on MS with Smatrix coefficients (NEW).\n"	;


	if(!parsecmdline("src_hyperspectral_fusion_L1_weights", arg, argc, argv, options, parameters))
		return EXIT_FAILURE;


	// Parameters
	int nlWeightsOpt = atoi(oNLwOpt.value);
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
		printf("ERROR :: src_hyperspectral_fusion_L1_weights :: uspampled HS image must be defined on high-resolution domain.\n");
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
		printf("ERROR :: src_hyperspectral_fusion_L1_weights :: spectral-downsampling matrix must be a one-channel image.\n");
		return EXIT_FAILURE;
	}

	if(Simg.w() != hs_channels)
	{
		printf("ERROR :: src_hyperspectral_fusion_L1_weights :: spectral-downsampling matrix width must be the number of HS channels.\n");
		return EXIT_FAILURE;
	}

	if(Simg.h() != ms_channels)
	{
		printf("ERROR :: src_hyperspectral_fusion_L1_weights :: spectral-downsampling matrix height must be the number of MS channels.\n");
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

	if(nlWeightsOpt == 1) // Weights computed on MS data (the same for all bands)
	{
		std::vector< std::vector <float> > wxy;
		std::vector< std::vector <int> > posxy;
		std::vector< std::vector <float> > wyx;
		std::vector< std::vector <int> > posyx;
		std::vector< std::vector <int> > posw;

		libUSTGHYPER::nlweights_ms(fM, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind, compwind, flagNormalize,
				ms_channels, width,  height);

		libUSTG::flimage outWxy(numNeighbours, dim);
		libUSTG::flimage outWyx(numNeighbours, dim);
		libUSTG::flimage outPosxy(numNeighbours, dim);
		libUSTG::flimage outPosyx(numNeighbours, dim);
		libUSTG::flimage outPosw(numNeighbours, dim);

		for(int i = 0; i < dim; i++)
		{
			int sizexy = (int) wxy[i].size();

			numNeighxy[i] = sizexy;

			for(int w = 0; w < sizexy; w++)
			{
				int index = i * numNeighbours + w;
				outWxy.v()[index] = wxy[i][w];
				outPosxy.v()[index] = posxy[i][w];
			}

			int sizeyx = (int) wyx[i].size();
			numNeighyx[i] = sizeyx;

			for(int w = 0; w < sizeyx; w++)
			{
				int index = i * numNeighbours + w;
				outWyx.v()[index] = wyx[i][w];
				outPosyx.v()[index] = posyx[i][w];
				outPosw.v()[index] = posw[i][w];
			}
		}

		outWxy.save(pWxy.value);
		outWyx.save(pWyx.value);
		outPosxy.save(pPosxy.value);
		outPosyx.save(pPosyx.value);
		outPosw.save(pPosw.value);

	} else if(nlWeightsOpt == 2) // Weights computed on MS data but balanced by the Smatrix coefficients
	{
		std::vector< std::vector< std::vector <float> > > wxy;
		std::vector< std::vector< std::vector <int> > > posxy;
		std::vector< std::vector< std::vector <float> > > wyx;
		std::vector< std::vector< std::vector <int> > > posyx;
		std::vector< std::vector< std::vector <int> > > posw;


		libUSTGHYPER:: nlweights_ms_Smatrix(fM, S, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind, compwind,
				flagNormalize, hs_channels, ms_channels, width, height);

		libUSTG::cflimage outWxy(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outWyx(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outPosxy(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outPosyx(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outPosw(numNeighbours, dim, hs_channels);

		for(int i = 0; i < dim; i++)
		{
			numNeighxy[i] = (int) wxy[0][i].size();
			numNeighyx[i] = (int) wyx[0][i].size();
		}

		for(int h = 0; h < hs_channels; h++)
		{
			for(int i = 0; i < dim; i++)
			{
				for(int w = 0; w < numNeighxy[i]; w++)
				{
					int index = i * numNeighbours + w;
					outWxy.v(h)[index] = wxy[h][i][w];
					outPosxy.v(h)[index] = posxy[h][i][w];
				}

				for(int w = 0; w < numNeighyx[i]; w++)
				{
					int index = i * numNeighbours + w;
					outWyx.v(h)[index] = wyx[h][i][w];
					outPosyx.v(h)[index] = posyx[h][i][w];
					outPosw.v(h)[index] = posw[h][i][w];
				}
			}
		}

		outWxy.save(pWxy.value);
		outWyx.save(pWyx.value);
		outPosxy.save(pPosxy.value);
		outPosyx.save(pPosyx.value);
		outPosw.save(pPosw.value);

	} else if(nlWeightsOpt == 3) // Weights computed on each upsampled HS band
	{
		std::vector< std::vector< std::vector <float> > > wxy;
		std::vector< std::vector< std::vector <int> > > posxy;
		std::vector< std::vector< std::vector <float> > > wyx;
		std::vector< std::vector< std::vector <int> > > posyx;
		std::vector< std::vector< std::vector <int> > > posw;

		libUSTGHYPER::nlweights_hs_band(fHint, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind, compwind,
				flagNormalize, hs_channels, width, height);

		libUSTG::cflimage outWxy(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outWyx(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outPosxy(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outPosyx(numNeighbours, dim, hs_channels);
		libUSTG::cflimage outPosw(numNeighbours, dim, hs_channels);

		for(int i = 0; i < dim; i++)
		{
			numNeighxy[i] = (int) wxy[0][i].size();
			numNeighyx[i] = (int) wyx[0][i].size();
		}

		for(int h = 0; h < hs_channels; h++)
		{
			for(int i = 0; i < dim; i++)
			{
				for(int w = 0; w < numNeighxy[i]; w++)
				{
					int index = i * numNeighbours + w;
					outWxy.v(h)[index] = wxy[h][i][w];
					outPosxy.v(h)[index] = posxy[h][i][w];
				}

				for(int w = 0; w < numNeighyx[i]; w++)
				{
					int index = i * numNeighbours + w;
					outWyx.v(h)[index] = wyx[h][i][w];
					outPosyx.v(h)[index] = posyx[h][i][w];
					outPosw.v(h)[index] = posw[h][i][w];
				}
			}
		}

		outWxy.save(pWxy.value);
		outWyx.save(pWyx.value);
		outPosxy.save(pPosxy.value);
		outPosyx.save(pPosyx.value);
		outPosw.save(pPosw.value);

	} else if(nlWeightsOpt == 4) // Weights computed on full upsampled HS data (the same for all bands)
	{
		std::vector< std::vector <float> > wxy;
		std::vector< std::vector <int> > posxy;
		std::vector< std::vector <float> > wyx;
		std::vector< std::vector <int> > posyx;
		std::vector< std::vector <int> > posw;

		libUSTGHYPER::nlweights_hs_3Dall(himgInt, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind,
				compwind, flagNormalize);

		libUSTG::flimage outWxy(numNeighbours, dim);
		libUSTG::flimage outWyx(numNeighbours, dim);
		libUSTG::flimage outPosxy(numNeighbours, dim);
		libUSTG::flimage outPosyx(numNeighbours, dim);
		libUSTG::flimage outPosw(numNeighbours, dim);

		for(int i = 0; i < dim; i++)
		{
			int sizexy = (int) wxy[i].size();
			numNeighxy[i] = sizexy;

			for(int w = 0; w < sizexy; w++)
			{
				int index = i * numNeighbours + w;
				outWxy.v()[index] = wxy[i][w];
				outPosxy.v()[index] = posxy[i][w];
			}

			int sizeyx = (int) wyx[i].size();
			numNeighyx[i] = sizeyx;

			for(int w = 0; w < sizeyx; w++)
			{
				int index = i * numNeighbours + w;
				outWyx.v()[index] = wyx[i][w];
				outPosyx.v()[index] = posyx[i][w];
				outPosw.v()[index] = posw[i][w];
			}
		}

		outWxy.save(pWxy.value);
		outWyx.save(pWyx.value);
		outPosxy.save(pPosxy.value);
		outPosyx.save(pPosyx.value);
		outPosw.save(pPosw.value);

	} else if(nlWeightsOpt == 5) // Weights computed on nearest upsampled HS bands
	{
		std::vector< std::vector <float> > wxy;
		std::vector< std::vector <int> > posxy;
		std::vector< std::vector <float> > wyx;
		std::vector< std::vector <int> > posyx;
		std::vector< std::vector <int> > posw;

		libUSTGHYPER::nl_weights_hs_3Dnearest(himgInt, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind,
				compwind, bandwind, flagNormalize);

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

	} else if(nlWeightsOpt == 6) // Neighbors selected on full upsampled HS data and weights computed just on associated HS band
	{
		std::vector< std::vector <float> > wxy;
		std::vector< std::vector <int> > posxy;
		std::vector< std::vector <float> > wyx;
		std::vector< std::vector <int> > posyx;
		std::vector< std::vector <int> > posw;

		libUSTGHYPER::nlweights_hs_3Dclassify_2Dband(himgInt, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours,
				reswind, compwind, flagNormalize);

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

	} else if(nlWeightsOpt == 7) // Neighbors selected on full MS data and weights computed just on associated HS band
	{
		std::vector< std::vector <float> > wxy;
		std::vector< std::vector <int> > posxy;
		std::vector< std::vector <float> > wyx;
		std::vector< std::vector <int> > posyx;
		std::vector< std::vector <int> > posw;

		libUSTGHYPER::nlweights_ms_3Dclassify_2Dband(himgInt, mimg, wxy, posxy, wyx, posyx, posw, hPatch,
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

	} else if(nlWeightsOpt == 8) // Neighbors selected on full MS data with Smatrix and weights computed just on associated HS band
	{
		std::vector< std::vector <float> > wxy;
		std::vector< std::vector <int> > posxy;
		std::vector< std::vector <float> > wyx;
		std::vector< std::vector <int> > posyx;
		std::vector< std::vector <int> > posw;

		libUSTGHYPER::nlweights_ms_3Dclassify_Smatrix_2Dband(himgInt, mimg, S, wxy, posxy, wyx, posyx, posw, hPatch,
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

	}else if(nlWeightsOpt == 9) // Neighbors selected on full MS data with Smatrix and weights computed just on associated HS band (new code)
	{
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

	}else
	{
		printf("ERROR :: src_hyperspectral_fusion :: incorrect option for computing NL weights:\n"
				"      - 1 : weights computed on MS (the same for all bands).\n"
				"      - 2 : weights computed on MS with Smatrix coefficients.\n"
				"      - 3 : weights computed on each upsampled HS band.\n"
				"      - 4 : weights computed on full upsampled HS (the same for all bands).\n"
				"      - 5 : weights computed on nearest upsampled HS bands.\n"
				"      - 6 : neighbors selected on full upsampled HS and weights computed just on associated HS band.\n"
				"      - 7 : neighbors selected on full MS and weights computed just on associated HS band.\n"
				"      - 8 : neighbors selected on full MS with Smatrix coefficients and weights computed just on associated HS band.");
		return EXIT_FAILURE;
	}


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
