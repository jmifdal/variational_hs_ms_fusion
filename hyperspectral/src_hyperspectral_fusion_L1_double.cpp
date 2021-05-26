#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "../library/libImage.h"
#include "../library/libHyperspectralDOUBLE.h"


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

    /*ParStruct pWxy_L2 = {"wxy", NULL, "input wxy NL weights L2 norm"};
    parameters.push_back(&pWxy_L2);
    ParStruct pWyx_L2 = {"wyx", NULL, "input wxy NL weights L2 norm"};
    parameters.push_back(&pWyx_L2);
    ParStruct pPosxy_L2 = {"posxy", NULL, "input posxy NL weights L2 norm"};
    parameters.push_back(&pPosxy_L2);
    ParStruct pPosyx_L2 = {"posyx", NULL, "input posyx NL weights L2 norm"};
    parameters.push_back(&pPosyx_L2);
    ParStruct pPosw_L2 = {"posw", NULL, "input posw NL weights L2 norm"};
    parameters.push_back(&pPosw_L2);
    ParStruct pNumNeighxy_L2 = {"numNeighbours", NULL, "input number of neighbours per pixel in wxy weights L2 norm"};
    parameters.push_back(&pNumNeighxy_L2);
    ParStruct pNumNeighyx_L2 = {"numNeighbours", NULL, "input number of neighbours per pixel in wyx weights L2 norm"};
    parameters.push_back(&pNumNeighyx_L2);*/


    char* arg = "INFO: Fusion of a high-resolution multispectral image and a low-resolution hyperspectral image by NL regularization,\n"
                "      L2 data-fidelity terms and L1 data-fitting constraint.\n"
                "      The algorithm only works for BLABLA Gaussian deconvolution and the standard deviation is required.\n";

    if(!parsecmdline("src_hyperspectral_fusion_L1", arg, argc, argv, options, parameters))
        return EXIT_FAILURE;


    // Parameters
    double lmbM = atof(oLmbM.value);
    double lmbH = atof(oLmbH.value);
    double mu = atof(oMu.value);
    double tol = atof(oTol.value);
    int maxIter = atoi(oMaxIter.value);
    double tau = atof(oTau.value);
    double sigma = atof(oSigma.value);
    double stdBlur = atof(pStdBlur.value);
    int sampling = atoi(pSampling.value);


    // MS and HS images
    libUSTGDOUBLE::cflimageDOUBLE mimg, himg;
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
        printf("ERROR :: src_hyperspectral_fusion_L1 :: MS and HS sizes must be related according to sampling factor.\n");
        return EXIT_FAILURE;
    }

    //Copying the multispectral and the hyperspectral data in a 2d matrix
    //the number of rows is the number of spectral channels, and each row contains
    //the pixels of the corresponding channel
    double **fM = new double*[ms_channels];
    for(int m = 0; m < ms_channels; m++)
    {
        fM[m] = new double[dim];
        libUSTGDOUBLE::fpCopy(mimg.v(m), fM[m], dim);
    }

    double **fH = new double*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        fH[h] = new double[s_dim];
        libUSTGDOUBLE::fpCopy(himg.v(h), fH[h], s_dim);
    }


    // Upsampled HS image
    libUSTGDOUBLE::cflimageDOUBLE himgInt;
    himgInt.load(pHyperInt.value);

    if((himgInt.w() != width) || (himgInt.h() != height))
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: uspampled HS image must be defined on high-resolution domain.\n");
        return EXIT_FAILURE;
    }

    if(himgInt.c() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: upsampled HS image must have as channels as HS data.\n");
        return EXIT_FAILURE;
    }

    double **fHint = new double*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        fHint[h] = new double[dim];
        libUSTGDOUBLE::fpCopy(himgInt.v(h), fHint[h], dim);
    }


    // Panchromatic images
    libUSTGDOUBLE::cflimageDOUBLE pimg, pimgInt;
    pimg.load(pHyperPan.value);
    pimgInt.load(pHyperPanInt.value);

    if((pimg.w() != width) || (pimg.h() != height) || (pimgInt.w() != width) || (pimgInt.h() != height))
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: panchromatic images must be defined on high-resolution domain.\n");
        return EXIT_FAILURE;
    }

    if((pimg.c() != hs_channels) || (pimgInt.c() != hs_channels))
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: panchromatic images must have as channels as HS data.\n");
        return EXIT_FAILURE;
    }

    //copying the panchromatic and the panchromatic from interpolated hyperspectral image
    //in a 2d structure
    double **PH = new double*[hs_channels];
    double **PHint = new double*[hs_channels];
    for(int h = 0; h < hs_channels; h++)
    {
        PH[h] = new double[dim];
        libUSTGDOUBLE::fpCopy(pimg.v(h), PH[h], dim);

        PHint[h] = new double[dim];
        libUSTGDOUBLE::fpCopy(pimgInt.v(h), PHint[h], dim);
    }


    // Spectral-downsampling operator
    libUSTGDOUBLE::flimageDOUBLE Simg;
    Simg.load(pSmatrix.value);

    if(Simg.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: spectral-downsampling matrix must be a one-channel image.\n");
        return EXIT_FAILURE;
    }

    if(Simg.w() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: spectral-downsampling matrix width must be the number of HS channels.\n");
        return EXIT_FAILURE;
    }

    if(Simg.h() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: spectral-downsampling matrix height must be the number of MS channels.\n");
        return EXIT_FAILURE;
    }

    //Copying S and S_transposed in a 2d structure
    double **S = new double*[ms_channels];

    for(int m = 0; m < ms_channels; m++)
    {
        S[m] = new double[hs_channels];

        for(int h = 0; h < hs_channels; h++)
            S[m][h] = Simg[m * hs_channels + h];
    }

    // Spectral-upsampling operator
    libUSTGDOUBLE::flimageDOUBLE Stimg;
    Stimg.load(pSmatrixT.value);

    if(Stimg.c() != 1)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: spectral-upsampling matrix must be a one-channel image.\n");
        return EXIT_FAILURE;
    }

    if(Stimg.w() != ms_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: spectral-upsampling matrix width must be the number of MS channels.\n");
        return EXIT_FAILURE;
    }

    if(Stimg.h() != hs_channels)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: spectral-upsampling matrix height must be the number of HS channels.\n");
        return EXIT_FAILURE;
    }

    double **St = new double*[hs_channels];

    for(int h = 0; h < hs_channels; h++)
    {
        St[h] = new double[ms_channels];

        for(int m = 0; m < ms_channels; m++)
            St[h][m] = Stimg[h * ms_channels + m];
    }

    // NL weights and apply fusion algorithm
    double **u = new double*[hs_channels];

    for(int h = 0; h < hs_channels; h++)
        u[h] = new double[dim];

    libUSTGDOUBLE::flimageDOUBLE imnumxy, imnumyx;
    imnumxy.load(pNumNeighxy.value);
    imnumyx.load(pNumNeighyx.value);

    int *numNeighxy = new int[dim];
    int *numNeighyx = new int[dim];

    //DEBUGG
    /*libUSTG::flimage imnumxy_L2, imnumyx_L2;
    imnumxy_L2.load(pNumNeighxy_L2.value);
    imnumyx_L2.load(pNumNeighyx_L2.value);

    int *numNeighxy_L2 = new int[dim];
    int *numNeighyx_L2 = new int[dim];*/



    for(int i = 0; i < dim; i++)
    {
        numNeighxy[i] = (int) imnumxy.v()[i];
        numNeighyx[i] = (int) imnumyx.v()[i];

        //DEBUGG
        /*numNeighxy_L2[i] = (int) imnumxy_L2.v()[i];
        numNeighyx_L2[i] = (int) imnumyx_L2.v()[i];*/
    }

    libUSTGDOUBLE::cflimageDOUBLE imwxy, imposxy, imwyx, imposyx, imposw;
    imwxy.load(pWxy.value);
    imposxy.load(pPosxy.value);
    imwyx.load(pWyx.value);
    imposyx.load(pPosyx.value);
    imposw.load(pPosw.value);

    //DEBUGG
    /*libUSTG::cflimage  imwxy_L2,imposxy_L2,imwyx_L2,imposyx_L2,imposw_L2;
    imwxy_L2.load(pWxy_L2.value);
    imposxy_L2.load(pPosxy_L2.value);
    imwyx_L2.load(pWyx_L2.value);
    imposyx_L2.load(pPosyx_L2.value);
    imposw_L2.load(pPosw_L2.value);*/


    int NLneighbours = imwxy.w();
    int NLdim = imwxy.h();
    int NLchannels = imwxy.c();

    //DEBUGG
    /*int NLneighbours_L2 = imwxy_L2.w();
    int NLdim_L2 = imwxy_L2.h();
    int NLchannels_L2 = imwxy_L2.c();*/

    if(NLdim != dim)
    {
        printf("ERROR :: src_hyperspectral_fusion_L1 :: the width in weights must be the number of high-resolution pixels.\n");
        return EXIT_FAILURE;
    }

    if(NLchannels == 1)
    {
        /*std::vector< std::vector <double> > wxy;
        std::vector< std::vector <int> > posxy;
        std::vector< std::vector <double> > wyx;
        std::vector< std::vector <int> > posyx;
        std::vector< std::vector <int> > posw;

        for(int i = 0; i < dim; i++)
        {
            std::vector <double> wxy_i;
            std::vector <int> posxy_i;
            std::vector <double> wyx_i;
            std::vector <int> posyx_i;
            std::vector <int> posw_i;

            for(int w = 0; w < numNeighxy[i]; w++)
            {
                int index = i * NLneighbours + w;
                wxy_i.push_back((double) imwxy.v()[index]);
                posxy_i.push_back((int) imposxy.v()[index]);
            }

            for(int w = 0; w < numNeighyx[i]; w++)
            {
                int index = i * NLneighbours + w;
                wyx_i.push_back((double) imwyx.v()[index]);
                posyx_i.push_back((int) imposyx.v()[index]);
                posw_i.push_back((int) imposw.v()[index]);
            }

            wxy.push_back(wxy_i);
            posxy.push_back(posxy_i);
            wyx.push_back(wyx_i);
            posyx.push_back(posyx_i);
            posw.push_back(posw_i);
        }

        libUSTGHYPERDOUBLE::hyperfusion_NL_L1(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau,
                                        sigma, tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);*/

    	printf("(CONVERSION TO DOUBLE) rien à faire \n");

    } else
    {
        if(NLchannels != hs_channels)
        {
            printf("ERROR :: src_hyperspectral_fusion_L1 :: the number of channels in weights must be the number of HS bands.\n");
            return EXIT_FAILURE;

        } else
        {
            std::vector< std::vector< std::vector <double> > > wxy;
            std::vector< std::vector< std::vector <int> > > posxy;
            std::vector< std::vector< std::vector <double> > > wyx;
            std::vector< std::vector< std::vector <int> > > posyx;
            std::vector< std::vector< std::vector <int> > > posw;

            //DEBUGG
            /*std::vector< std::vector< std::vector <float> > > wxy_L2;
            std::vector< std::vector< std::vector <int> > > posxy_L2;
            std::vector< std::vector< std::vector <float> > > wyx_L2;
            std::vector< std::vector< std::vector <int> > > posyx_L2;
            std::vector< std::vector< std::vector <int> > > posw_L2;*/

            for(int h = 0; h < hs_channels; h++)
            {
                std::vector< std::vector <double> > wxy_h;
                std::vector< std::vector <int> > posxy_h;
                std::vector< std::vector <double> > wyx_h;
                std::vector< std::vector <int> > posyx_h;
                std::vector< std::vector <int> > posw_h;

                //DEBUGG
                /*std::vector< std::vector <float> > wxy_h_L2;
                std::vector< std::vector <int> > posxy_h_L2;
                std::vector< std::vector <float> > wyx_h_L2;
                std::vector< std::vector <int> > posyx_h_L2;
                std::vector< std::vector <int> > posw_h_L2;*/

                for(int i = 0; i < dim; i++)
                {
                    std::vector <double> wxy_i;
                    std::vector <int> posxy_i;
                    std::vector <double> wyx_i;
                    std::vector <int> posyx_i;
                    std::vector <int> posw_i;

                    //DEBUGG
                    /*std::vector <float> wxy_i_L2;
                    std::vector <int> posxy_i_L2;
                    std::vector <float> wyx_i_L2;
                    std::vector <int> posyx_i_L2;
                    std::vector <int> posw_i_L2;*/



                    for(int w = 0; w < numNeighxy[i]; w++)
                    {
                        int index = i * NLneighbours + w;
                        wxy_i.push_back((double) imwxy.v(h)[index]);
                        posxy_i.push_back((int) imposxy.v(h)[index]);

                        //DEBUGG
                        /*int index_L2 = i * NLneighbours_L2 + w;
                        wxy_i_L2.push_back((float) imwxy_L2.v(h)[index_L2]);
                        posxy_i_L2.push_back((int) imposxy_L2.v(h)[index_L2]);*/
                    }

                    for(int w = 0; w < numNeighyx[i]; w++)
                    {
                        int index = i * NLneighbours + w;
                        wyx_i.push_back((double) imwyx.v(h)[index]);
                        posyx_i.push_back((int) imposyx.v(h)[index]);
                        posw_i.push_back((int) imposw.v(h)[index]);

                        //DEBUGG
                        /*int index_L2 = i * NLneighbours_L2 + w;
                        wyx_i_L2.push_back((float) imwyx_L2.v(h)[index_L2]);
                        posyx_i_L2.push_back((int) imposyx_L2.v(h)[index_L2]);
                        posw_i_L2.push_back((int) imposw_L2.v(h)[index_L2]);*/
                    }

                    wxy_h.push_back(wxy_i);
                    posxy_h.push_back(posxy_i);
                    wyx_h.push_back(wyx_i);
                    posyx_h.push_back(posyx_i);
                    posw_h.push_back(posw_i);

                    //DEBUGG
                    /*wxy_h_L2.push_back(wxy_i_L2);
                    posxy_h_L2.push_back(posxy_i_L2);
                    wyx_h_L2.push_back(wyx_i_L2);
                    posyx_h_L2.push_back(posyx_i_L2);
                    posw_h_L2.push_back(posw_i_L2);*/
                }

                wxy.push_back(wxy_h);
                posxy.push_back(posxy_h);
                wyx.push_back(wyx_h);
                posyx.push_back(posyx_h);
                posw.push_back(posw_h);

                //DEBUGG
                /*wxy_L2.push_back(wxy_h_L2);
                posxy_L2.push_back(posxy_h_L2);
                wyx_L2.push_back(wyx_h_L2);
                posyx_L2.push_back(posyx_h_L2);
                posw_L2.push_back(posw_h_L2);*/
            }

            //printf("Hello \n");

            libUSTGHYPERDOUBLE::hyperfusion_NL_L1_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau,
                                                           sigma, tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height);

            //DEBUGG
            /*libUSTGHYPER::hyperfusion_NL_L1_wh(u, fH, fM, S, St, PH, PHint, fHint, wxy, posxy, wyx, posyx, posw, lmbH, lmbM, mu, tau,
                                               sigma, tol, maxIter, sampling, stdBlur, hs_channels, ms_channels, width, height,
                                               wxy_L2, posxy_L2, wyx_L2, posyx_L2, posw_L2);*/
            //printf("Hello again \n");
        }
    }


    // Save output
    libUSTGDOUBLE::cflimageDOUBLE output(width, height, hs_channels);

    for(int h = 0; h < hs_channels; h++)
    	libUSTGDOUBLE::fpCopy(u[h], output.v(h), dim);

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
