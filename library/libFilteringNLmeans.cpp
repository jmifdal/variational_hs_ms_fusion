#include "libFilteringNLmeans.h"



namespace libUSTG
{
	
	
    
    
    
	void fpNeighFiltering(int iDWin,			// Half size of comparison window
                          int iDBloc,		// Half size of research window
                          float fFiltSpa,	// Noise parameter
                          float fFiltIm,     // Filtering parameter
                          float **fpI,		// Input
                          float **fpA,		// Auxiliary input for distance computation
                          float **fpO,		// Output
                          float *imask,      // Mask of pixels to be filtered
                          int iType,
                          int iIChannels, int iAChannels, int iWidth,int iHeight)
	{
		
		
		int iwxh = iWidth * iHeight;                //! length of each channel
	 	
		//! exponential tabulation
		int iLutLength = (int) rintf((float) LUTMAX * (float) LUTPRECISION);
		float *fpLut = new float[iLutLength];
		wxFillExpLut(fpLut,iLutLength);
        
        
        //! filtering parameters
		float fImH2 = fFiltIm * fFiltIm;
		float fSpaH2 = fFiltSpa * fFiltSpa;
		
        
		//! clear output
		for (int ii=0; ii < iIChannels; ii++) fpCopy(fpI[ii], fpO[ii], iwxh);
		
		
		// PROCESS STARTS
		// for each pixel (x,y)
//#pragma omp parallel shared(fpI, fpA, fpO, fpLut)
		{
			
			
//#pragma omp for schedule(dynamic) nowait
			
			for(int y=0; y < iHeight ; y++)
			{
				
				//! denoised value at a certain pixel
				float *fpDenoised = new float[iIChannels];
				
                float **fpV = new float*[iIChannels];
                float **fpW = new float*[iIChannels];
                for (int ii=0; ii < iIChannels; ii++)
                {
                    fpV[ii] = new float[ (2*iDBloc+1) * (2*iDBloc+1)];
                    fpW[ii] = new float[ (2*iDBloc+1) * (2*iDBloc+1)];
                }
                
                
				for(int x=0 ; x < iWidth; x++)
                    if ( imask[ y * iWidth + x] > 0.0f )
                    {
                        
                        int ixy = y*iWidth+x;
                        
                        //! reduce the size of the comparison window if we are near the boundary
                        int iDWin0 = MIN(iDWin,MIN(iWidth-1-x,MIN(iHeight-1-y,MIN(x,y))));
                        
                        //! research zone depending on the boundary and the size of the window
                        int imin=MAX(x-iDBloc,iDWin0);
                        int jmin=MAX(y-iDBloc,iDWin0);
                        
                        int imax=MIN(x+iDBloc,iWidth-1-iDWin0);
                        int jmax=MIN(y+iDBloc,iHeight-1-iDWin0);
                        
                        
                        //!  clear current denoised value
                        fpClear(fpDenoised, 0.0f, iIChannels);
                        
                        
                        float fMaxWeight = 0.0f;            	//! maximum of weights for central pixel
                        float fTotalWeight = 0.0f;      		//! sum of weights
                        float fWeight = 0.0f;
                        
                        
                        //!  compute l2 distances
                        int iV = 0;
                        for(int j=jmin; j <= jmax; j++)
                            for(int i=imin ; i <= imax; i++, iV++)
                                if (i!=x || j!=y)
                                {
                                    
                                    int iij = j*iWidth+i;
                                    
                                    float fDif = 0.0f;
                                    if (fpA)
                                        fDif = fiL2FloatDist(fpA,fpA,x,y,i,j,iDWin0, iDWin0, iAChannels,iWidth,iWidth);
                                    else
                                        fDif = fiL2FloatDist(fpI,fpI,x,y,i,j,iDWin0, iDWin0,  iIChannels,iWidth,iWidth);
                                    
                                    fDif = fDif / fImH2;
                                    fWeight = wxSLUT(fDif,fpLut);
                                    
                                    float fSpa;
                                    if (fFiltSpa > 0.0f)
                                    {
                                        fSpa= (float) (i-x) * (float) (i-x) + (float) (j-y) * (float) (j-y);
                                        fSpa = fSpa / fSpaH2;
                                        fWeight *= wxSLUT(fSpa,fpLut);
                                    }
                                    
                                    
                                    if (fWeight > fMaxWeight) fMaxWeight = fWeight;
                                    fTotalWeight += fWeight;
                                    
                                    
                                    for (int ii=0; ii < iIChannels; ii++)  fpDenoised[ii] +=  fWeight * fpI[ii][iij];
                                    
                                    for (int ii=0; ii < iIChannels; ii++)
                                    {
                                        fpV[ii][iV] =  fpI[ii][iij];
                                        fpW[ii][iV] =  fWeight;
                                    }
                                    
                                    
                                }
                        
                        
                        
                        //! Accumulate Central Point
                        for (int ii=0; ii < iIChannels; ii++)  fpDenoised[ii] +=  fMaxWeight * fpI[ii][ixy];
                        fTotalWeight += fMaxWeight;
                        
                        for (int ii=0; ii < iIChannels; ii++)
                        {
                            fpV[ii][iV] =  fpI[ii][ixy];
                            fpW[ii][iV] =  fMaxWeight;
                        }
                        iV++;
                        
                        
                        if (fTotalWeight > fTiny)
                        {
                            if (iType == 0)
                            {
                                for (int ii=0; ii < iIChannels; ii++)  fpO[ii][y*iWidth+x] = fpDenoised[ii] / fTotalWeight;
							}
                            else if (iType==1)
                            {
                                
                                float hfTotalWeight = 0.5 * fTotalWeight;
                                
                                for (int ii=0; ii < iIChannels; ii++)
                                {
                                    
                                    fpQuickSort(fpV[ii], fpW[ii], iV);
                                    
                                    float pweight = 0.0f;
                                    int jj = 0;
                                    for (jj=0; jj < iV && pweight < hfTotalWeight; jj++)
                                    {
                                        pweight += fpW[ii][jj];
                                        
                                    }
                                    
                                    fpO[ii][y*iWidth+x] = fpV[ii][jj];
                                    
                                }
                                
                            }
                            else if (iType==2)
                            {
                            
                                for (int ii=0; ii < iIChannels; ii++)
                                {
                                    float fMax= -fLarge;
                                    float fValue = 0.0;
                                    for (int jj=0; jj < iV; jj++)
                                    {
                                    
                                        if ( fpW[ii][jj] * fpV[ii][jj] > fMax) {fMax = fpW[ii][jj] * fpV[ii][jj];fValue=fpV[ii][jj];}
                                    
                                    }
                                    
                                    fpO[ii][y*iWidth+x] = fValue;

                                }
                                
                            } else
                            {
                                
                                for (int ii=0; ii < iIChannels; ii++)
                                {
                                    float fMin= fLarge;
                                    float fValue = 0.0;
                                    for (int jj=0; jj < iV; jj++)
                                    {
                                        
                                        if ( fpW[ii][jj] * fpV[ii][jj] < fMin) {fMin = fpW[ii][jj] * fpV[ii][jj];fValue=fpV[ii][jj];}
                                        
                                    }
                                    
                                    fpO[ii][y*iWidth+x] = fValue;
                                    
                                }

                                
                            }
                            
                            
                        }
                        else
                        {
                            
                            for (int ii=0; ii < iIChannels; ii++)  fpO[ii][y*iWidth+x] =  fpI[ii][y*iWidth+x];
                        }
                        
                        
                    }
				
                
                delete[] fpDenoised;
                for (int ii=0; ii < iIChannels; ii++)
                {
                    delete[] fpV[ii];
                    delete[] fpW[ii];
                }
                
                delete[] fpV;
                delete[] fpW;
                
			}
			
			
		}
		
		
		
		
	}
	
    
    
    
    
    
    
    void compute_weights(int iDWin,         // Half size of comparison window
                         int iDBloc,		// Half size of research window
                         float fFiltSpa,	// Noise parameter
                         float fFiltIm,     // Filtering parameter
                         float **fpI,		// Input
                         float *imask,      // Mask of pixels to be filtered
                         float **fpWeight,
                         int iIChannels,  int iWidth,int iHeight)
    {
        
    	//! exponential tabulation
		int iLutLength = (int) rintf((float) LUTMAX * (float) LUTPRECISION);
		float *fpLut = new float[iLutLength];
		wxFillExpLut(fpLut,iLutLength);
        
        
        //! filtering parameters
		float fImH2 = fFiltIm * fFiltIm;
		float fSpaH2 = fFiltSpa * fFiltSpa;
        
        
//#pragma omp parallel shared(fpI, fpLut, fpWeight)
		{
			
			
//#pragma omp for schedule(dynamic) nowait
			
			for(int y=0; y < iHeight ; y++)
			{
				
                
				for(int x=0 ; x < iWidth; x++)
                    if ( imask[ y * iWidth + x] > 0.0f )
                    {
                        
                        int ixy = y*iWidth+x;
                        
                        fpWeight[ixy] = new float[(2*iDBloc+1) * (2*iDBloc+1)];
                        
                        //! reduce the size of the comparison window if we are near the boundary
                        int iDWin0 = MIN(iDWin,MIN(iWidth-1-x,MIN(iHeight-1-y,MIN(x,y))));
                        
                        //! research zone depending on the boundary and the size of the window
                        int imin=MAX(x-iDBloc,iDWin0);
                        int jmin=MAX(y-iDBloc,iDWin0);
                        
                        int imax=MIN(x+iDBloc,iWidth-1-iDWin0);
                        int jmax=MIN(y+iDBloc,iHeight-1-iDWin0);
                        
                        float fMaxWeight = 0.0f;            	//! maximum of weights for central pixel
                        float fWeight = 0.0f;
                        
                        
                        //!  compute l2 distances
                        int iV=0;
                        int iCentral = 0;
                        for(int j=jmin; j <= jmax; j++)
                            for(int i=imin ; i <= imax; i++, iV++)
                                if (i!=x || j!=y)
                                {
                                    
                                    float fDif = 0.0f;
                                    fDif = fiL2FloatDist(fpI,fpI,x,y,i,j,iDWin0, iDWin0,  iIChannels,iWidth,iWidth);
                                    
                                    fDif = fDif / fImH2;
                                    fWeight = wxSLUT(fDif,fpLut);
                                    
                                    float fSpa;
                                    if (fFiltSpa > 0.0f)
                                    {
                                        fSpa= (float) (i-x) * (float) (i-x) + (float) (j-y) * (float) (j-y);
                                        fSpa = fSpa / fSpaH2;
                                        fWeight *= wxSLUT(fSpa,fpLut);
                                    }
                                    
                                    
                                    if (fWeight > fMaxWeight) fMaxWeight = fWeight;
                                    
                                    fpWeight[ixy][iV] = fWeight;
                                    
                                } else iCentral = iV;
                        
                        fpWeight[ixy][iCentral] = fMaxWeight;
                    }
                
                
            }
            
        }
        
        
    }
    
    
    
    void fpNeighFilteringIt(int iDWin,		// Half size of comparison window
                            int iDBloc,		// Half size of research window
                            float fFiltSpa,	// Noise parameter
                            float fFiltIm,    // Filtering parameter
                            int   nIter,        // Number of iterations
                            float **fpI,		// Input
                            float **fpA,		// Auxiliary input for distance computation
                            float **fpO,		// Output
                            float *imask,     // Mask of pixels to be filtered
                            int iType,
                            int iIChannels, int iAChannels, int iWidth,int iHeight)
	{
		
		
		int iwxh = iWidth * iHeight;                //! length of each channel
	 	float **fpWeights = new float*[iwxh];
        
        
        
        if (fpA)
            compute_weights(iDWin, iDBloc, fFiltSpa, fFiltIm, fpA, imask,fpWeights, iAChannels, iWidth,iHeight);
        else
            compute_weights(iDWin, iDBloc, fFiltSpa, fFiltIm, fpI, imask,fpWeights, iIChannels, iWidth,iHeight);
        
		
        
		//! clear output
		for (int ii=0; ii < iIChannels; ii++) fpCopy(fpI[ii], fpO[ii], iwxh);
        
        for (int it = 0; it < nIter; it++)
        {
            
            
            // PROCESS STARTS
            // for each pixel (x,y)
//#pragma omp parallel shared(fpI, fpA, fpO, fpWeights)
            {
                
                
//#pragma omp for schedule(dynamic) nowait
                
                for(int y=0; y < iHeight ; y++)
                {
                    
                    //! denoised value at a certain pixel
                    float *fpDenoised = new float[iIChannels];
                    
                    float **fpV = new float*[iIChannels];
                    float **fpW = new float*[iIChannels];
                    for (int ii=0; ii < iIChannels; ii++)
                    {
                        fpV[ii] = new float[ (2*iDBloc+1) * (2*iDBloc+1)];
                        fpW[ii] = new float[ (2*iDBloc+1) * (2*iDBloc+1)];
                    }
                    
                    
                    for(int x=0 ; x < iWidth; x++)
                        if ( imask[ y * iWidth + x] > 0.0f )
                        {
                            
                            int ixy = y*iWidth+x;
                            
                            //! reduce the size of the comparison window if we are near the boundary
                            int iDWin0 = MIN(iDWin,MIN(iWidth-1-x,MIN(iHeight-1-y,MIN(x,y))));
                            
                            //! research zone depending on the boundary and the size of the window
                            int imin=MAX(x-iDBloc,iDWin0);
                            int jmin=MAX(y-iDBloc,iDWin0);
                            
                            int imax=MIN(x+iDBloc,iWidth-1-iDWin0);
                            int jmax=MIN(y+iDBloc,iHeight-1-iDWin0);
                            
                            
                            //!  clear current denoised value
                            fpClear(fpDenoised, 0.0f, iIChannels);
                            
                            
                         
                            float fTotalWeight = 0.0f;      		//! sum of weights
                            float fWeight = 0.0f;
                            
                            
                            //!  compute l2 distances
                            int iV = 0;
                            for(int j=jmin; j <= jmax; j++)
                                for(int i=imin ; i <= imax; i++, iV++)
                                {
                                    
                                    int iij = j*iWidth+i;
                                    
                                    fWeight = fpWeights[ixy][iV];
                                    
                                    for (int ii=0; ii < iIChannels; ii++)  fpDenoised[ii] +=  fWeight * fpI[ii][iij];
                                    
                                    fTotalWeight += fWeight;
                           
                                    
                                    for (int ii=0; ii < iIChannels; ii++)
                                    {
                                        fpV[ii][iV] =  fpI[ii][iij];
                                        fpW[ii][iV] =  fWeight;
                                    }
                                    
                                    
                                }
                            
                            
                        
                            
                            if (fTotalWeight > fTiny)
                            {
                                if (iType == 0)
                                {
                                    for (int ii=0; ii < iIChannels; ii++)  fpO[ii][y*iWidth+x] = fpDenoised[ii] / fTotalWeight;
                                    
                                }else if (iType==1)
                                {
                                    
                                    float hfTotalWeight = 0.5 * fTotalWeight;
                                    
                                    for (int ii=0; ii < iIChannels; ii++)
                                    {
                                        
                                        fpQuickSort(fpV[ii], fpW[ii], iV);
                                        
                                        float pweight = 0.0f;
                                        int jj = 0;
                                        for (jj=0; jj < iV && pweight < hfTotalWeight; jj++)
                                        {
                                            pweight += fpW[ii][jj];
                                            
                                        }
                                        
                                        fpO[ii][y*iWidth+x] = fpV[ii][jj];
                                        
                                    }
                                    
                                }
                                else if (iType==2)
                                {
                                    
                                    for (int ii=0; ii < iIChannels; ii++)
                                    {
                                        float fMax= -fLarge;
                                        float fValue = 0.0f;
                                        
                                        for (int jj=0; jj < iV; jj++)
                                        {
                                            
                                            if ( fpW[ii][jj] * fpV[ii][jj] > fMax) {fMax = fpW[ii][jj] * fpV[ii][jj]; fValue=fpV[ii][jj];}
                                            
                                        }
                                        
                                        fpO[ii][y*iWidth+x] = fValue;
                                        
                                    }
                                    
                                } else
                                {
                                    
                                    for (int ii=0; ii < iIChannels; ii++)
                                    {
                                        float fMin= fLarge;
                                        float fValue = 0.0f;
                                        for (int jj=0; jj < iV; jj++)
                                        {
                                            
                                            if ( fpW[ii][jj] * fpV[ii][jj] < fMin) {fMin = fpW[ii][jj] * fpV[ii][jj];fValue=fpV[ii][jj];}
                                            
                                        }
                                        
                                        fpO[ii][y*iWidth+x] = fValue;
                                        
                                    }
                                    
                                    
                                }
                                
                                
                                
                                
                            }
                            
                        }
                    
                    
                    delete[] fpDenoised;
                    for (int ii=0; ii < iIChannels; ii++)
                    {
                        delete[] fpV[ii];
                        delete[] fpW[ii];
                    }
                    
                    delete[] fpV;
                    delete[] fpW;
                    
                }
                
                
            }
            
            for (int ii=0; ii < iIChannels; ii++) fpCopy(fpO[ii], fpI[ii], iwxh);
            
            
		}
        
        
		
		
	}
    
    
    
    
    
    
    
}

