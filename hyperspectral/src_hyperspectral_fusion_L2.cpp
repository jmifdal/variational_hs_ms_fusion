#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "../library/libImage.h"
#include "../library/libHyperspectral.h"


int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
    OptStruct oLmbM = {"m:", 0, "1", NULL, "trade-off parameter for multispectral data term"}; options.push_back(&oLmbM);
    OptStruct oLmbH = {"h:", 0, "1", NULL, "trade-off parameter for hyperspectral data term"}; options.push_back(&oLmbH);
    OptStruct oMu = {"u:", 0, "1", NULL, "trade-off paramater for radiometric constraint"}; options.push_back(&oMu);
    OptStruct oTol = {"e:", 0, "0.000001", NULL, "stopping precision of primal-dual algorithm "}; options.push_back(&oTol);
    OptStruct oMaxIter = {"i:", 0, "1000", NULL, "max allowed iterations for primal-dual algorithm"}; options.push_back(&oMaxIter);
    OptStruct oTau = {"t:", 0, "0.05", NULL, "step-size of primal-dual gradient descent algorithm"}; options.push_back(&oTau);
    OptStruct oSigma = {"s:", 0, "0.05", NULL, "step-size of primal-dual gradient ascent algorithm"}; options.push_back(&oSigma);

    std::vector<ParStruct *> parameters;
    ParStruct pMulti = {"multispectral", NULL, "input high-resolution multispectral image"};
    parameters.push_back(&pMulti);
    ParStruct pHyper = {"hyperspectral", NULL, "input low-resolution hyperspectral image"};
    parameters.push_back(&pHyper);
    ParStruct pHyperInt = {"hyperspectral_int", NULL, "input spatially upsampled hyperspectral image"};
    parameters.push_back(&pHyperInt);
    ParStruct pHyperPan = {"hyper_pan", NULL, "input pan computed from multispectral bands using Smatrix"};
    parameters.push_back(&pHyperPan);
    ParStruct pHyperPanInt = {"hyper_pan_int", NULL, "input pan computed from upsampled lowres multispectral bands using Smatrix"};
    parameters.push_back(&pHyperPanInt);
    ParStruct pWxy = {"wxy", NULL, "input wxy NL weights"};
    parameters.push_back(&pWxy);
    ParStruct pWyx = {"wyx", NULL, "input wxy NL weights"};
    parameters.push_back(&pWyx);
    ParStruct pPosxy = {"posxy", NULL, "input posxy NL weights"};
    parameters.push_back(&pPosxy);
    ParStruct pPosyx = {"posyx", NULL, "input posyx NL weights"};
    parameters.push_back(&pPosyx);
    ParStruct pPosw = {"posw", NULL, "input posw NL weights"};
    parameters.push_back(&pPosw);
    ParStruct pNumNeighxy = {"numNeighbours", NULL, "input number of neighbours per pixel in wxy weights"};
    parameters.push_back(&pNumNeighxy);
    ParStruct pNumNeighyx = {"numNeighbours", NULL, "input number of neighbours per pixel in wyx weights"};
    parameters.push_back(&pNumNeighyx);
    ParStruct pSmatrix = {"Smatrix", NULL, "input spectral-downsampling matrix operator"};
    parameters.push_back(&pSmatrix);
    ParStruct pSmatrixT = {"SmatrixT", NULL, "input spectral-upsampling matrix operator (transposed to spectral-downsampling matrix)"};
    parameters.push_back(&pSmatrixT);
    ParStruct pStdBlur = {"stdBlur", NULL, "input standard deviation of Gaussian kernel for spatial downsampling"};
    parameters.push_back(&pStdBlur);
    ParStruct pSampling = {"sampling", NULL, "input sampling factor for spatial downsampling"};
    parameters.push_back(&pSampling);
    ParStruct pOutput = {"output", NULL, "output high-resolution hyperspectral image"};
    parameters.push_back(&pOutput);

    char* arg = "INFO: Fusion of a high-resolution multispectral image and a low-resolution hyperspectral image by NL regularization,\n"
                "      L2 data-fidelity terms and L2 data-fitting constraint.\n"
                "      The algorithm only works for Gaussian deconvolution and the standard deviation is required.\n";

    if(!parsecmdline("src_hyperspectral_fusion_L2", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;


    // Parameters
    float lmbM = atof(oLmbM.value);
    float lmbH = atof(oLmbH.value);
    float mu = atof(oMu.value);
    float tol = atof(oTol.value);
    int maxIter = atoi(oMaxIter.value);
    float tau = atof(oTau.value);
    float sigma = atof(oSigma.value);
    float stdBlur = atof(pStdBlur.value);
    int sampling = atoi(pSampling.value);


    // MS and HS images
    libUSTG::cflimage mimg, himg;
    mimg.load(pMulti.value);
    himg.load(pHyper.value);

    int width = mimg.w();
    int height = mimg.h();
    int ms_channels = mimg.c();
    int dim = width * height;

    int s_width = himg.w();
    int s_height = himg.h();
    int hs_channels = himg.c();
    int s_dim = s_width * s_height;

    if((width / sampling != s_width) || (height / sampling != s_height))
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: MS and HS sizes must be related according to sampling factor.\n");
        return EXIT_FAILURE;
    }

    //Copying the multispectral and the hyperspectral data in a 2d matrix
    //the number of rows is the number of spectral channels, and each row contains
    //the pixels of the corresponding channel
    float **fM = new float*[ms_channels];
    for(int m = 0; m < ms_channels; m++)
    {
        fM[m] = new float[dim];
        libUSTG::fpCopy(mimg.v(m), fM[m], dim);
    }

    float **fH = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        fH[h] = new float[s_dim];
        libUSTG::fpCopy(himg.v(h), fH[h], s_dim);
    }


    // Upsampled HS image
    libUSTG::cflimage himgInt;
    himgInt.load(pHyperInt.value);

    if((himgInt.w() != width) || (himgInt.h() != height))
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: uspampled HS image must be defined on high-resolution domain.\n");
        return EXIT_FAILURE;
    }

    if(himgInt.c() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: upsampled HS image must have as channels as HS data.\n");
        return EXIT_FAILURE;
    }

    float **fHint = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        fHint[h] = new float[dim];
        libUSTG::fpCopy(himgInt.v(h), fHint[h], dim);
    }


    // Panchromatic images
    libUSTG::cflimage pimg, pimgInt;
    pimg.load(pHyperPan.value);
    pimgInt.load(pHyperPanInt.value);

    if((pimg.w() != width) || (pimg.h() != height) || (pimgInt.w() != width) || (pimgInt.h() != height))
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: panchromatic images must be defined on high-resolution domain.\n");
        return EXIT_FAILURE;
    }

    if((pimg.c() != hs_channels) || (pimgInt.c() != hs_channels))
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: panchromatic images must have as channels as HS data.\n");
        return EXIT_FAILURE;
    }

    //copying the panchromatic and the panchromatic from interpolated hyperspectral image
    //in a 2d structure
    float **PH = new float*[hs_channels];
    float **PHint = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        PH[h] = new float[dim];
        libUSTG::fpCopy(pimg.v(h), PH[h], dim);

        PHint[h] = new float[dim];
        libUSTG::fpCopy(pimgInt.v(h), PHint[h], dim);
    }


    // Spectral-downsampling operator
    libUSTG::flimage Simg;
    Simg.load(pSmatrix.value);

    if(Simg.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: spectral-downsampling matrix must be a one-channel image.\n");
        return EXIT_FAILURE;
    }

    if(Simg.w() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: spectral-downsampling matrix width must be the number of HS channels.\n");
        return EXIT_FAILURE;
    }

    if(Simg.h() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: spectral-downsampling matrix height must be the number of MS channels.\n");
        return EXIT_FAILURE;
    }

    //Copying S and S_transposed in a 2d structure
    float **S = new float*[ms_channels];

    for(int m = 0; m < ms_channels; m++)
    {
        S[m] = new float[hs_channels];

        for(int h = 0; h < hs_channels; h++)
            S[m][h] = Simg[m * hs_channels + h];
    }

    // Spectral-upsampling operator
    libUSTG::flimage Stimg;
    Stimg.load(pSmatrixT.value);

    if(Stimg.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: spectral-upsampling matrix must be a one-channel image.\n");
        return EXIT_FAILURE;
    }

    if(Stimg.w() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: spectral-upsampling matrix width must be the number of MS channels.\n");
        return EXIT_FAILURE;
    }

    if(Stimg.h() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: spectral-upsampling matrix height must be the number of HS channels.\n");
        return EXIT_FAILURE;
    }

    float **St = new float*[hs_channels];

    for(int h = 0; h < hs_channels; h++)
    {
        St[h] = new float[ms_channels];

        for(int m = 0; m < ms_channels; m++)
            St[h][m] = Stimg[h * ms_channels + m];
    }

    // NL weights and apply fusion algorithm
    float **u = new float*[hs_channels];

    for(int h = 0; h < hs_channels; h++)
        u[h] = new float[dim];

    libUSTG::flimage imnumxy, imnumyx;
    imnumxy.load(pNumNeighxy.value);
    imnumyx.load(pNumNeighyx.value);

    int *numNeighxy = new int[dim];
    int *numNeighyx = new int[dim];

    for(int i = 0; i < dim; i++)
    {
        numNeighxy[i] = (int) imnumxy.v()[i];
        numNeighyx[i] = (int) imnumyx.v()[i];
    }

    libUSTG::cflimage imwxy, imposxy, imwyx, imposyx, imposw;
    imwxy.load(pWxy.value);
    imposxy.load(pPosxy.value);
    imwyx.load(pWyx.value);
    imposyx.load(pPosyx.value);
    imposw.load(pPosw.value);

    int NLneighbours = imwxy.w();
    int NLdim = imwxy.h();
    int NLchannels = imwxy.c();

    if(NLdim != dim)
    {
        printf("ERROR :: src_hyperspectral_fusion_L2 :: the width in weights must be the number of high-resolution pixels.\n");
        return EXIT_FAILURE;
    }

    if(NLchannels == 1)
    {
        std::vector< std::vector <float> > wxy;
        std::vector< std::vector <int> > posxy;
        std::vector< std::vector <float> > wyx;
        std::vector< std::vector <int> > posyx;
        std::vector< std::vector <int> > posw;

        for(int i = 0; i < dim; i++)
        {
            std::vector <float> wxy_i;
            std::vector <int> posxy_i;
            std::vector <float> wyx_i;
            std::vector <int> posyx_i;
            std::vector <int> posw_i;

            for(int w = 0; w < numNeighxy[i]; w++)
            {
                int index = i * NLneighbours + w;
                wxy_i.push_back((float) imwxy.v()[index]);
                posxy_i.push_back((int) imposxy.v()[index]);
            }

            for(int w = 0; w < numNeighyx[i]; w++)
            {
                int index = i * NLneighbours + w;
                wyx_i.push_back((float) imwyx.v()[index]);
                posyx_i.push_back((int) imposyx.v()[index]);
                posw_i.push_back((int) imposw.v()[index]);
            }

            wxy.push_back(wxy_i);
            posxy.push_back(posxy_i);
            wyx.push_back(wyx_i);
            posyx.push_back(posyx_i);
            posw.push_back(posw_i);
        }

        libUSTGHYPER::hyperfusion_NL(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau,
                                        sigma, tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);

    } else
    {
        if(NLchannels != hs_channels)
        {
            printf("ERROR :: src_hyperspectral_fusion_L2 :: the number of channels in weights must be the number of HS bands.\n");
            return EXIT_FAILURE;

        } else
        {
            std::vector< std::vector< std::vector <float> > > wxy;
            std::vector< std::vector< std::vector <int> > > posxy;
            std::vector< std::vector< std::vector <float> > > wyx;
            std::vector< std::vector< std::vector <int> > > posyx;
            std::vector< std::vector< std::vector <int> > > posw;

            for(int h = 0; h < hs_channels; h++)
            {
                std::vector< std::vector <float> > wxy_h;
                std::vector< std::vector <int> > posxy_h;
                std::vector< std::vector <float> > wyx_h;
                std::vector< std::vector <int> > posyx_h;
                std::vector< std::vector <int> > posw_h;

                for(int i = 0; i < dim; i++)
                {
                    std::vector <float> wxy_i;
                    std::vector <int> posxy_i;
                    std::vector <float> wyx_i;
                    std::vector <int> posyx_i;
                    std::vector <int> posw_i;

                    for(int w = 0; w < numNeighxy[i]; w++)
                    {
                        int index = i * NLneighbours + w;
                        wxy_i.push_back((float) imwxy.v(h)[index]);
                        posxy_i.push_back((int) imposxy.v(h)[index]);
                    }

                    for(int w = 0; w < numNeighyx[i]; w++)
                    {
                        int index = i * NLneighbours + w;
                        wyx_i.push_back((float) imwyx.v(h)[index]);
                        posyx_i.push_back((int) imposyx.v(h)[index]);
                        posw_i.push_back((int) imposw.v(h)[index]);
                    }

                    wxy_h.push_back(wxy_i);
                    posxy_h.push_back(posxy_i);
                    wyx_h.push_back(wyx_i);
                    posyx_h.push_back(posyx_i);
                    posw_h.push_back(posw_i);
                }

                wxy.push_back(wxy_h);
                posxy.push_back(posxy_h);
                wyx.push_back(wyx_h);
                posyx.push_back(posyx_h);
                posw.push_back(posw_h);
            }

            //printf("Hello \n");

            libUSTGHYPER::hyperfusion_NL_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau,
                                               sigma, tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);
            //printf("Hello again \n");
        }
    }


    // Save output
    libUSTG::cflimage output(width, height, hs_channels);

    for(int h = 0; h < hs_channels; h++)
        libUSTG::fpCopy(u[h], output.v(h), dim);

    output.save(pOutput.value);


    // Delete allocated memory
    for(int m = 0; m < ms_channels; m++)
    {
        delete[] fM[m];
        delete[] S[m];
    }

    for(int h = 0; h < hs_channels; h++)
    {
        delete[] fH[h];
        delete[] fHint[h];
        delete[] PH[h];
        delete[] PHint[h];
        delete[] St[h];
        delete[] u[h];
    }

    delete[] fM;
    delete[] fH;
    delete[] fHint;
    delete[] PH;
    delete[] PHint;
    delete[] S;
    delete[] St;
    delete[] u;

    return EXIT_SUCCESS;
}
















/*************************************************
 * Previous code, saved on the 01/06/2019
 *************************************************/

/*#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "../library/libImage.h"
#include "../library/libHyperspectral.h"
#include "../library/libHyperspectralWeights.h"
//#include "./libHyperspectralWeights.h"


int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
    OptStruct oLmbM = {"m:", 0, "1", NULL, "trade-off parameter for multispectral data term"}; options.push_back(&oLmbM);
    OptStruct oLmbH = {"h:", 0, "1", NULL, "trade-off parameter for hyperspectral data term"}; options.push_back(&oLmbH);
    OptStruct oMu = {"u:", 0, "1", NULL, "trade-off paramater for radiometric constraint"}; options.push_back(&oMu);
    OptStruct oNLwOpt = {"o:", 0, "3", NULL, "NL weights on upsampled hyperspectral (1), multispectral (2) or multispectral * Smatrix (3)\n"}; options.push_back(&oNLwOpt);
    OptStruct oHPatchDist = {"c:", 0, "2.5", NULL, "filtering parameter for patch distance in NL weights"}; options.push_back(&oHPatchDist);
    OptStruct oHSpatDist = {"d:", 0, "2.5", NULL, "filtering parameter for spatial distance in NL weights"}; options.push_back(&oHSpatDist);
    OptStruct oNneig = {"n:", 0, "25", NULL, "number of (most similar) neighboring pixels used in NL weights"}; options.push_back(&oNneig);
    OptStruct oBandWindSize = {"b:", 0, "1", NULL, "half-size of spectral window in NL weights"}; options.push_back(&oBandWindSize);
    OptStruct oResWindSize = {"r:", 0, "7", NULL, "half-size of research window in NL weights"}; options.push_back(&oResWindSize);
    OptStruct oPatchSize = {"p:", 0, "1", NULL, "half-size of patches in NL weights"}; options.push_back(&oPatchSize);
    OptStruct oTol = {"e:", 0, "0.000001", NULL, "stopping precision of primal-dual algorithm "}; options.push_back(&oTol);
    OptStruct oMaxIter = {"i:", 0, "1000", NULL, "max allowed iterations for primal-dual algorithm"}; options.push_back(&oMaxIter);
    OptStruct oTau = {"t:", 0, "0.05", NULL, "step-size of primal-dual gradient descent algorithm"}; options.push_back(&oTau);
    OptStruct oSigma = {"s:", 0, "0.05", NULL, "step-size of primal-dual gradient ascent algorithm"}; options.push_back(&oSigma);
    
    std::vector<ParStruct *> parameters;
    ParStruct pMulti = {"multispectral", NULL, "input high-resolution multispectral image"}; parameters.push_back(&pMulti);
    ParStruct pHyper = {"hyperspectral", NULL, "input low-resolution hyperspectral image"}; parameters.push_back(&pHyper);
    ParStruct pHyperInt = {"hyperspectral_int", NULL, "input spatially upsampled hyperspectral image"}; parameters.push_back(&pHyperInt);
    ParStruct pHyperPan = {"hyper_pan", NULL, "input pan computed from multispectral bands using Smatrix"}; parameters.push_back(&pHyperPan);
    ParStruct pHyperPanInt = {"hyper_pan_int", NULL, "input pan computed from upsampled lowres multispectral bands using Smatrix"}; parameters.push_back(&pHyperPanInt);
    ParStruct pSmatrix = {"Smatrix", NULL, "input spectral-downsampling matrix operator"}; parameters.push_back(&pSmatrix);
    ParStruct pSmatrixT = {"SmatrixT", NULL, "input spectral-upsampling matrix operator (transposed to spectral-downsampling matrix)"}; parameters.push_back(&pSmatrixT);
    ParStruct pStdBlur = {"stdBlur", NULL, "input standard deviation of Gaussian kernel for spatial downsampling"}; parameters.push_back(&pStdBlur);
    ParStruct pSampling = {"sampling", NULL, "input sampling factor for spatial downsampling"}; parameters.push_back(&pSampling);
    ParStruct pOutput = {"output", NULL, "output high-resolution hyperspectral image"}; parameters.push_back(&pOutput);
   
    char* arg = "INFO: Fusion of a high-resolution multispectral image and a low-resolution hyperspectral image by NL regularization and L2 data terms.\n"
                "      The algorithm only works for Gaussian deconvolution and the standard deviation is required.\n"
                "      Nonlocal weights are computed according to the following options:\n"
                "        - 1 : weights computed on MS data (the same for all bands).\n"
                "        - 2 : weights computed on MS data but balanced by the Smatrix coefficients.\n"
                "        - 3 : weights computed on each upsampled HS band.\n"
                "        - 4 : weights computed on full upsampled HS data (the same for all bands).\n"
                "        - 5 : weights computed on nearest upsampled HS bands.\n"
                "        - 6 : neighbors selected on full upsampled HS data and weights computed just on associated HS band.\n";
    
    if(!parsecmdline("src_hyperspectral_fusion", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;
    
    
    // Parameters
    float lmbM = atof(oLmbM.value);
    float lmbH = atof(oLmbH.value);
    float mu = atof(oMu.value);
    int nlWeightsOpt = atoi(oNLwOpt.value);
    float hPatch = atof(oHPatchDist.value);
    float hSpatial = atof(oHSpatDist.value);
    int numNeighbours = atoi(oNneig.value);
    int bandwind = atoi(oBandWindSize.value);
    int reswind = atoi(oResWindSize.value);
    int compwind = atoi(oPatchSize.value);
    float tol = atof(oTol.value);
    int maxIter = atoi(oMaxIter.value);
    float tau = atof(oTau.value);
    float sigma = atof(oSigma.value);
    float stdBlur = atof(pStdBlur.value);
    int sampling = atoi(pSampling.value);
    
    // Image data
    libUSTG::cflimage mimg, himg;
    mimg.load(pMulti.value);
    himg.load(pHyper.value);
    
    int width = mimg.w();
    int height = mimg.h();
    int ms_channels = mimg.c();
    int dim = width * height;
    
    int s_width = himg.w();
    int s_height = himg.h();
    int hs_channels = himg.c();
    int s_dim = s_width * s_height;
    
    if((width / sampling != s_width) || (height / sampling != s_height))
    {
        printf("ERROR :: src_hyperspectral_fusion :: multispectral and hyperspectral sizes must be related according to sampling factor.\n");
        return EXIT_FAILURE;
    }
    
    float **fM = new float*[ms_channels];
    for(int m = 0; m < ms_channels; m++)
    {
        fM[m] = new float[dim];
        libUSTG::fpCopy(mimg.v(m), fM[m], dim);
    }
    
    float **fH = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        fH[h] = new float[s_dim];
        libUSTG::fpCopy(himg.v(h), fH[h], s_dim);
    }
    
    // Upsampled hyperspectral image
    libUSTG::cflimage himgInt;
    himgInt.load(pHyperInt.value);
    
    if((himgInt.w() != width) || (himgInt.h() != height))
    {
        printf("ERROR :: src_hyperspectral_fusion :: uspampled hyperspectral image must be defined on high-resolution domain.\n");
        return EXIT_FAILURE;
    }
    
    if(himgInt.c() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion :: upsampled hyperspectral image must have as channels as hyperspectral data.\n");
        return EXIT_FAILURE;
    }
    
    float **fHint = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        fHint[h] = new float[dim];
        libUSTG::fpCopy(himgInt.v(h), fHint[h], dim);
    }
    
    // Panchromatic images
    libUSTG::cflimage pimg, pimgInt;
    pimg.load(pHyperPan.value);
    pimgInt.load(pHyperPanInt.value);
    
    if((pimg.w() != width) || (pimg.h() != height) || (pimgInt.w() != width) || (pimgInt.h() != height))
    {
        printf("ERROR :: src_hyperspectral_fusion :: panchromatic images must be defined on high-resolution domain.\n");
        return EXIT_FAILURE;
    }
    
    if((pimg.c() != hs_channels) || (pimgInt.c() != hs_channels))
    {
        printf("ERROR :: src_hyperspectral_fusion :: panchromatic images must have as channels as hyperspectral data.\n");
        return EXIT_FAILURE;
    }
    
    float **PH = new float*[hs_channels];
    float **PHint = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        PH[h] = new float[dim];
        libUSTG::fpCopy(pimg.v(h), PH[h], dim);
        
        PHint[h] = new float[dim];
        libUSTG::fpCopy(pimgInt.v(h), PHint[h], dim);
    }
    
    // Spectral-downsampling operator
    libUSTG::flimage Simg;
    Simg.load(pSmatrix.value);
    
    if(Simg.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_fusion :: spectral-downsampling matrix must be a one-channel image.\n");
        return EXIT_FAILURE;
    }
    
    if(Simg.w() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion :: spectral-downsampling matrix width must be the number of hyperspectral channels.\n");
        return EXIT_FAILURE;
    }
    
    if(Simg.h() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion :: spectral-downsampling matrix height must be the number of multispectral channels.\n");
        return EXIT_FAILURE;
    }
    
    float **S = new float*[ms_channels];
    
    for(int m = 0; m < ms_channels; m++)
    {
        S[m] = new float[hs_channels];
        
        for(int h = 0; h < hs_channels; h++)
            S[m][h] = Simg[m * hs_channels + h];
    }
    
    // Spectral-upsampling operator
    libUSTG::flimage Stimg;
    Stimg.load(pSmatrixT.value);
    
    if(Stimg.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_fusion :: spectral-upsampling matrix must be a one-channel image.\n");
        return EXIT_FAILURE;
    }
    
    if(Stimg.w() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion:: spectral-upsampling matrix width must be the number of multispectral channels.\n");
        return EXIT_FAILURE;
    }
    
    if(Stimg.h() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion :: spectral-upsampling matrix height must be the number of hyperspectral channels.\n");
        return EXIT_FAILURE;
    }
    
    float **St = new float*[hs_channels];
    
    for(int h = 0; h < hs_channels; h++)
    {
        St[h] = new float[ms_channels];
        
        for(int m = 0; m < ms_channels; m++)
            St[h][m] = Stimg[h * ms_channels + m];
    }
    
    // Compute NL weights and apply fusion model
    float **u = new float*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
        u[h] = new float[dim];
    
    int flagNormalize = 1;
    
    if(nlWeightsOpt == 1) // Weights computed on MS data (the same for all bands)
    {
        // Compute NL weights
        std::vector< std::vector <float> > wxy;
        std::vector< std::vector <int> > posxy;
        std::vector< std::vector <float> > wyx;
        std::vector< std::vector <int> > posyx;
        std::vector< std::vector <int> > posw;
        
        libUSTGHYPER::nlweights_ms(fM, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind, compwind, flagNormalize,
                                   ms_channels, width,  height);
        
        // Fusion
        libUSTGHYPER::hyperfusion_NL(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau, sigma, tol,
                                     maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);
        
    } else if(nlWeightsOpt == 2) // Weights computed on MS data but balanced by the Smatrix coefficients
    {
        // Compute NL weights
        std::vector< std::vector< std::vector <float> > > wxy;
        std::vector< std::vector< std::vector <int> > > posxy;
        std::vector< std::vector< std::vector <float> > > wyx;
        std::vector< std::vector< std::vector <int> > > posyx;
        std::vector< std::vector< std::vector <int> > > posw;
        
        clock_t t1, t2;
        t1=t2=clock();
        
        libUSTGHYPER:: nlweights_ms_Smatrix(fM, S, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind, compwind,
                                            flagNormalize, hs_channels, ms_channels, width, height);
        
        t1=clock()-t1;
        // Fusion
        libUSTGHYPER::hyperfusion_NL_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau, sigma,
                                        tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);
        t2=clock()-t2;
        double time1 = ((double) t1) / CLOCKS_PER_SEC;
        double time2 = ((double) t2) / CLOCKS_PER_SEC;
        printf("%f %f\n", time1, time2);
    
    } else if(nlWeightsOpt == 3) // Weights computed on each upsampled HS band
    {
        // Compute NL weights
        std::vector< std::vector< std::vector <float> > > wxy;
        std::vector< std::vector< std::vector <int> > > posxy;
        std::vector< std::vector< std::vector <float> > > wyx;
        std::vector< std::vector< std::vector <int> > > posyx;
        std::vector< std::vector< std::vector <int> > > posw;
        
        libUSTGHYPER::nlweights_hs_band(fHint, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind, compwind,
                                        flagNormalize, hs_channels, width, height);
        
        // Fusion
        libUSTGHYPER::hyperfusion_NL_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau, sigma,
                                        tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);
        
    } else if(nlWeightsOpt == 4) // Weights computed on full upsampled HS data (the same for all bands)
    {
        // Compute NL weights
        std::vector< std::vector <float> > wxy;
        std::vector< std::vector <int> > posxy;
        std::vector< std::vector <float> > wyx;
        std::vector< std::vector <int> > posyx;
        std::vector< std::vector <int> > posw;
        
        libUSTGHYPER::nlweights_hs_3Dall(himgInt, wxy, posxy, wyx, posyx, posw, hPatch, hSpatial, numNeighbours, reswind,
                                         compwind, flagNormalize);
        
        // Fusion
        libUSTGHYPER::hyperfusion_NL(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau, sigma, tol,
                                     maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);
    
    } else if(nlWeightsOpt == 5) // Weights computed on nearest upsampled HS bands
    {
    	// Compute NL weights
    	std::vector< std::vector <float> > wxy0;
    	std::vector< std::vector <int> > posxy0;
    	std::vector< std::vector <float> > wyx0;
    	std::vector< std::vector <int> > posyx0;
    	std::vector< std::vector <int> > posw0;

    	libUSTGHYPER::nl_weights_hs_3Dnearest(himgInt, wxy0, posxy0, wyx0, posyx0, posw0, hPatch, hSpatial, numNeighbours, reswind, compwind,
    			bandwind, flagNormalize);

    	// Declaring arrays for weights and positions
    	std::vector< std::vector< std::vector <float> > > wxy;
    	std::vector< std::vector< std::vector <int> > > posxy;
    	std::vector< std::vector< std::vector <float> > > wyx;
    	std::vector< std::vector< std::vector <int> > > posyx;
    	std::vector< std::vector< std::vector <int> > > posw;

    	//conversion of the weights from a 2D to a 3D structure
    	//looping through the rows of the 2D structure
    	for(int k=0;k<=(hs_channels*dim - dim);k=k+dim)
    	{

    		std::vector< std::vector <float> > wxy_h;
    		std::vector< std::vector <int> > posxy_h;
    		std::vector< std::vector <float> > wyx_h;
    		std::vector< std::vector <int> > posyx_h;
    		std::vector< std::vector <int> > posw_h;

    		//looping through the pixels of each channel
    		for(int i=k;i<=k+dim-1;i++)
    		{

    			std::vector <float> wxy_i;
    			std::vector <int> posxy_i;
    			std::vector <float> wyx_i;
    			std::vector <int> posyx_i;
    			std::vector <int> posw_i;

    			//looping through the weights of the neighbors of each pixel i
    			for(int j=0;j<(int)wxy0[i].size();j++)
    			{
    				wxy_i.push_back(wxy0[i][j]);
    				posxy_i.push_back(posxy0[i][j]);

    			}

    			for(int j=0;j<(int)wyx0[i].size();j++)
    			{
    				wyx_i.push_back(wyx0[i][j]);
    				posyx_i.push_back(posyx0[i][j]);
    				posw_i.push_back(posw0[i][j]);
    			}

    			wxy_h.push_back(wxy_i);
    			posxy_h.push_back(posxy_i);
    			wyx_h.push_back(wyx_i);
    			posyx_h.push_back(posyx_i);
    			posw_h.push_back(posw_i);


    		}
    		wxy.push_back(wxy_h);
    		posxy.push_back(posxy_h);
    		wyx.push_back(wyx_h);
    		posyx.push_back(posyx_h);
    		posw.push_back(posw_h);

    	}

    	wxy0.clear(); posxy0.clear(); wyx0.clear(); posyx0.clear(); posw0.clear();

    	// Fusion
    	libUSTGHYPER::hyperfusion_NL_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau, sigma,
    			tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);

    } else if(nlWeightsOpt == 6) // Neighbors selected on full upsampled HS data and weights computed just on associated HS band
    {
        // Compute NL weights
        std::vector< std::vector <float> > wxy0;
        std::vector< std::vector <int> > posxy0;
        std::vector< std::vector <float> > wyx0;
        std::vector< std::vector <int> > posyx0;
        std::vector< std::vector <int> > posw0;
        
        libUSTGHYPER::nlweights_hs_3Dclassify_2Dband(himgInt, wxy0, posxy0, wyx0, posyx0, posw0, hPatch, hSpatial, numNeighbours, reswind,
                                                     compwind, flagNormalize);
        
        // Resize NL weights
        std::vector< std::vector< std::vector <float> > > wxy;
        std::vector< std::vector< std::vector <int> > > posxy;
        std::vector< std::vector< std::vector <float> > > wyx;
        std::vector< std::vector< std::vector <int> > > posyx;
        std::vector< std::vector< std::vector <int> > > posw;
        
        for(int h = 0; h < hs_channels; h++)
        {
            std::vector< std::vector <float> > wxy_h;
            std::vector< std::vector <int> > posxy_h;
            std::vector< std::vector <float> > wyx_h;
            std::vector< std::vector <int> > posyx_h;
            std::vector< std::vector <int> > posw_h;
            
            for(int i = 0; i < dim; i++)
            {
                std::vector <float> wxy_i;
                std::vector <int> posxy_i;
                std::vector <float> wyx_i;
                std::vector <int> posyx_i;
                std::vector <int> posw_i;
                
                for(int w = 0; w < (int) wxy0[h*dim+i].size(); w++)
                {
                    wxy_i.push_back(wxy0[h*dim+i][w]);
                    posxy_i.push_back(posxy0[h*dim+i][w]);
                    wyx_i.push_back(wyx0[h*dim+i][w]);
                    posyx_i.push_back(posyx0[h*dim+i][w]);
                    posw_i.push_back(posw0[h*dim+i][w]);
                }
                
                wxy_h.push_back(wxy_i);
                posxy_h.push_back(posxy_i);
                wyx_h.push_back(wyx_i);
                posyx_h.push_back(posyx_i);
                posw_h.push_back(posw_i);
            }
            
            wxy.push_back(wxy_h);
            posxy.push_back(posxy_h);
            wyx.push_back(wyx_h);
            posyx.push_back(posyx_h);
            posw.push_back(posw_h);
        }
        
        wxy0.clear(); posxy0.clear(); wyx0.clear(); posyx0.clear(); posw0.clear();
        
        // Fusion
        libUSTGHYPER::hyperfusion_NL_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau, sigma,
                                        tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);
        
    } else
    {
        printf("ERROR :: src_hyperspectral_fusion :: incorrect option for computing NL weights:\n"
               "        - 1 : weights computed on MS data (the same for all bands).\n"
               "        - 2 : weights computed on MS data but balanced by the Smatrix coefficients.\n"
               "        - 3 : weights computed on each upsampled HS band.\n"
               "        - 4 : weights computed on full upsampled HS data (the same for all bands).\n"
               "        - 5 : weights computed on nearest upsampled HS bands.\n"
               "        - 6 : neighbors selected on full upsampled HS data and weights computed just on associated HS band.");
        return EXIT_FAILURE;
    }


    // Save output
    libUSTG::cflimage output(width, height, hs_channels);
    
    for(int h = 0; h < hs_channels; h++)
        libUSTG::fpCopy(u[h], output.v(h), dim);
    
    output.save(pOutput.value);
    
    
    // Delete allocated memory
    for(int m = 0; m < ms_channels; m++)
    {
        delete[] fM[m];
        delete[] S[m];
    }
    
    for(int h = 0; h < hs_channels; h++)
    {
        delete[] fH[h];
        delete[] fHint[h];
        delete[] PH[h];
        delete[] PHint[h];
        delete[] St[h];
        delete[] u[h];
    }
    
    delete[] fM;
    delete[] fH;
    delete[] fHint;
    delete[] PH;
    delete[] PHint;
    delete[] S;
    delete[] St;
    delete[] u;
                                        
    return EXIT_SUCCESS;
}
*/
