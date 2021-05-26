#ifndef _LIBDENOISING_H_
#define _LIBDENOISING_H_

#include "libBasic.h"


#define MNPCA 2
#define VERBOSE 0



namespace libUSTG
{
	

    
    
	
    void fpNeighFiltering(int iDWin,			// Half size of comparison window
                          int iDBloc,           // Half size of research window
                          float fFiltSpa,       // Noise parameter
                          float fFiltIm,        // Filtering parameter
                          float **fpI,          // Input
                          float **fpA,          // Auxiliary input for distance computation
                          float **fpO,          // Output
                          float *imask,         // Mask of pixels to be filtered
                          int iType,            // iType=0 (average), 1 (median), 2 (max), 3 (min)
                          int iIChannels, int iAChannels, int iWidth,int iHeight);
    
    
    void fpNeighFilteringIt(int iDWin,          // Half size of comparison window
                            int iDBloc,         // Half size of research window
                            float fFiltSpa,     // Noise parameter
                            float fFiltIm,      // Filtering parameter
                            int   nIter,        // Number of iterations
                            float **fpI,		// Input
                            float **fpA,		// Auxiliary input for distance computation
                            float **fpO,		// Output
                            float *imask,       // Mask of pixels to be filtered
                            int iType,          // iType=0 (average), 1 (median), 2 (max), 3 (min)
                            int iIChannels, int iAChannels, int iWidth,int iHeight);
	
}


#endif
