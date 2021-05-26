#ifndef _library__
#define _library__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <time.h>
#include <cassert>
#include <vector>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <omp.h>
#include <fftw3.h>
#include <algorithm>

namespace libUSTG
{
    
    
#ifndef PI
#define PI 3.14159265358979323846264338327
#endif
    
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
    
#define dTiny 1e-10
#define fTiny 0.0000000001f
#define fLarge 10000000000.0f
#define dLarge 1e+10
    
#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )
    
    
    
    
    //! JL: For fast convolutions
    struct kernel_item {
        float v;
        int i, j;
    };
    
    struct sorted_kernel {
        int nitems;
        int w, h;
        int ncut;
        struct kernel_item *items;
    };
    
    
    
    //! Forward declarations
    class laMatrix;
    class laVector;
    
    
    
    
    void ustg_exit(const char* message);
    
    
    
    void ipClear(int *ipI, int iValue, int iLength);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Float Value operations
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    void fpClear(float *fpI,float fValue, int iLength);
	void fpCopy(float *fpI,float *fpO, int iLength);
	
	float fpMax(float *u,int *pos, int size);
	float fpMin(float *u,int *pos,int size);
	
    
	float fpVar(float *u,int size);
	float fpMean(float *u,int size);
    float fpMedian(float *u,int size);
    
    
    void fpCombine(float *u,float a,float *v,float b, float *w,  int size);
    
    
    void fiImageDrawCircle(float *igray, int pi,int pj, float radius, float value, int width, int height);
    void fiImageDrawLine(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height);
    
    
    void fpBinarize(float *u, float *v, float value, int inverse, int size);
    
    
    float fpDistLp(float *u, float *v, int p, int size);
    float fpDistLp(float *u, float *v, float *m, int p, int size);
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Double Value operations
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    double dpVar(double *u,int size);
	double dpMean(double *u,int size);
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Float pointer ordering
	//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
	void fpQuickSort(float *fpI, int iLength, int inverse = 0 );
	void fpQuickSort(float *fpI, float *fpO, int iLength, int inverse = 0);
	
	
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	//! Image Conversion
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
#define COEFF_YR 0.299
#define COEFF_YG 0.587
#define COEFF_YB 0.114
	
    void fiRgb2Yuv(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height);
	void fiYuv2Rgb(float *r,float *g,float *b,float *y,float *u,float *v, int width,int height);
	
	void fiRgb2YuvO(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height);
	void fiYuvO2Rgb(float *r,float *g,float *b,float *y,float *u,float *v, int width,int height);
    
    void fiRgb2Conrat(float *R, float *G, float *B, float *L, float *C1, float *C2, int dim);
    void fiConrat2Rgb(float *R, float *G, float *B, float *L, float *C1, float *C2, int dim);
    
    void fiRgb2LenzCarmona(float *R, float *G, float *B, float *L, float *C1, float *C2, int dim);
    void fiLenzCarmona2Rgb(float *R, float *G, float *B, float *L, float *C1, float *C2, int dim);
    
    void fiRgb2Ycbcr(float *R, float *G, float *B, float *Y, float *Cr, float *Cb, int dim);
    void fiYcbcr2Rgb(float *Y, float *Cb, float *Cr, float *R, float *G, float *B, int dim);
    
    void fiRgb2Hsv(float *R, float *G, float *B, float *H, float *S, float *V, int dim);
    void fiHsv2Rgb(float *H, float *S, float *V, float *R, float *G, float *B, int dim);
    
    void fiRgb2Hsl(float *R, float *G, float *B, float *H, float *S, float *L, int dim);
    void fiHsl2Rgb(float *H, float *S, float *L, float *R, float *G, float *B, int dim);
    
    void fiRgb2Hsi(float *R, float *G, float *B, float *H, float *S, float *I, int dim);
    void fiHsi2Rgb(float *H, float *S, float *I, float *R, float *G, float *B, int dim);
	
    //! XYZ color of the D65 white point */
#define WHITEPOINT_X	0.950456
#define WHITEPOINT_Y	1.0
#define WHITEPOINT_Z	1.088754
    
    //! sRGB gamma correction, transforms R to R'
#define GAMMACORRECTION(t)	\
(((t) <= 0.0031306684425005883) ? \
(12.92*(t)) : (1.055*pow((t), 0.416666666666666667) - 0.055))
    
    //! Inverse sRGB gamma correction, transforms R' to R
#define INVGAMMACORRECTION(t)	\
(((t) <= 0.0404482362771076) ? \
((t)/12.92) : pow(((t) + 0.055)/1.055, 2.4))
    
    //! CIE-L*a*b* f function (used to convert XYZ to L*a*b*)
#define LABF(t)	\
((t >= 8.85645167903563082e-3) ? \
pow(t,0.333333333333333) : (841.0/108.0)*(t) + (4.0/29.0))
    
    //! CIE-L*a*b* inverse f function
#define LABINVF(t)	\
((t >= 0.206896551724137931) ? \
((t)*(t)*(t)) : (108.0/841.0)*((t) - (4.0/29.0)))
    
    void fiRgb2Xyz(float *R, float *G, float *B, float *X, float *Y, float *Z, int dim);
    void fiXyz2Rgb(float *X, float *Y, float *Z, float *R, float *G, float *B, int dim);
    void fiXyz2Lab(float *X, float *Y, float *Z, float *L, float *a, float *b, int dim);
    void fiLab2Xyz(float *L, float *a, float *b, float *X, float *Y, float *Z, int dim);
    void fiRgb2Lab(float *R, float *G, float *B, float *L, float *a, float *b, int dim);
    void fiLab2Rgb(float *L, float *a, float *b, float *R, float *G, float *B, int dim);
    
	
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	//! Image Convolution
	//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    
#define BOUNDARY_CONDITION_NEUMANN 0
#define BOUNDARY_CONDITION_SYMMETRIC 1
    
	float* fiFloatGaussKernel(float std, int & size);
	void fiFloatDirectionalGaussKernel(float xsigma, float ysigma, float angle, float *kernel, int kwidth, int kheight);
    float * fiFloatDirectionalGaussKernelS(float xsigma, float ysigma, float angle, float *kernel, int kwidth, int kheight, int sign);

    
    void fiFloatHorizontalConvolution(float *u, float *v, int width, int height, float *kernel, int ksize, int boundary);
    void fiFloatVerticalConvolution(float *u, float *v, int width, int height, float *kernel,int ksize, int boundary);
	void fiGaussianConvol(float *u, float *v, int width, int height, float sigma, int boundary);
	void fiConvol(float *u,float *v,int width,int height,float *kernel,int kwidth,int kheight, int boundary);
	void fiSepConvol(float *u,float *v,int width,int height,float *xkernel, int xksize, float *ykernel, int yksize, int boundary);
	
    //! JL Lisani
    void sort_kernel(float *u, int w, int h, struct sorted_kernel *skernel, float pkernel);
    void fiConvol_skernel(float *u,float *v,int width,int height, struct sorted_kernel *skernel, int boundary);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	//! Tabulation
	//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
#define LUTMAX 30.0
#define dLUTMAX 3e+1
#define LUTMAXM1 29.0
#define dLUTMAXM1 29e+0
#define LUTPRECISION 1000.0
#define dLUTPRECISION 1e+3
	
	
	void  wxFillExpLut(float *lut, int size);
    void  wxFillExpLut(double *lut, int size);
	float wxSLUT(float dif, float *lut);
    double wxSLUT(double dif, double *lut);
    
    
    
    
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	//! Patch Distances
	//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
	float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int width0, int width1);
    float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int width0, int width1);
    float fiL2FloatDist(float **u0,float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int channels, int width0, int width1);
	
    
    float fiBCFloatDist(float **u0,float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int channels, int width0, int width1);
    float fiCFloatDist(float **u0,float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int channels, int width0, int width1);
    
    
	float fiL2FloatWDist ( float * u0, float *u1, int i0, int j0,int i1,int j1,int xradius, int yradius, float *kernel, int width0, int width1);
    float fiL2FloatWDist ( float ** u0, float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, float *kernel, int channels, int width0, int width1);
	
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //! Sampling functions
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void fiImageSample(float *igray,float *ogray, int factor, int width, int height);

    float *fiImageSample(float *input, float sampling_factor, int high_width,
                         int high_height, int & low_width, int & low_height);
    
    void fiImageSampleAglomeration(float *igray,float *ogray, int factor, int width, int height);
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //! Interpolation functions
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    float bicubic_interpolation_at(
                                   const float *input, // image to be interpolated
                                   const float  uu,    // x component of the vector field
                                   const float  vv,    // y component of the vector field
                                   const int    nx,    // image width
                                   const int    ny,    // image height
                                   const int    binterpolation,  // interpolate boundary if activated, else return bvalue
                                   const float  bvalue // return this value outside the region
    );
    
    
    void bicubic_interpolation_zoom(const float *input, const int nx, const int ny,  const float  fFactor, const int binterpolation, const float bvalue, float *output );
    void bicubic_interpolation_warp( const float *input, const float *u,  const float *v,  const int nx,  const int ny,  const int binterpolation, const float  bvalue, float *output);
    void nn_interpolation_zoom(const float *input, const int    nx, const int    ny,  const float  fFactor, float *output);
    
    void bicubic_homography_interpolation(float *input, int width, int height, laMatrix &H, float bg, float *out, int nwidth, int nheight);

    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //! FFT Stuff
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    void fft1d(float *Xr,float *Xi,float *Yr,float *Yi, int inverse, int n);
	void fft2d(float *in_re, float *in_im,float *out_re, float *out_im, int i_flag, int width, int height);
	void fiFFTShearPascal(float dAmountOfShear, float *pFloatImageInput,float *pFloatImageOutput, int iAxis, float fDelta, int PleaseDoNotCenter, int width, int height);
	void fiFFTZoom(float *in, float *out, float zoomx, float zoomy, int width, int height);
	
    
    //! Convolution
    void fft_Gaussian_convolution(float **convolved, float **data, float std, int num_channels, int width, int height);
    void fft_Gaussian_convolution(float *convolved, float *data, float std, int width, int height);
    void fft_convolution(float *convolved, float *data, float std, int width, int height);
    
    
    
 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
    //! Patch Statistics
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void fiComputeIntegralImage(float *in, float *out, int width, int height);
   
    void fiPatchStatistics(float *fpIn, float *fpMinV, float *fpMaxV, float *fpMeanV, float *fpVarV, float *fpMedianV, float fRadius, int iWidth, int iHeight);
    void fiPatchMin(float *fpIn, float *fpMinV, float fRadius, int iWidth, int iHeight);
    void fiPatchMax(float *fpIn, float *fpMaxV, float fRadius, int iWidth, int iHeight);
    void fiPatchMean(float *fpIn, float *fpMeanV, float fRadius, int iWidth, int iHeight);
    void fiPatchVar(float *fpIn, float *fpVarV, float fRadius, int iWidth, int iHeight);
    void fiPatchMedian(float *fpIn, float *fpMedianV, float fRadius, int iWidth, int iHeight);
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //! Splines Stuff  OLD
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    float vMW(float *in,int x,int y,float bg, int width, int height);
    void keysMW(float *c,float t,float a);
    void spline3MW(float *c,float t);
    void init_splinenMW(float *a,int n);
    float ipowMW(float x,int n);
    void splinenMW(float *c,float t,float *a,int n);
    float initcausalMW(float *c,int n,float z);
    float initanticausalMW(float *c,int n,float z);
    void invspline1DMW(float *c,int size,float *z,int npoles);
	void finvsplineMW(float *in,int order,float *out, int width, int height);
    float evaluate_splineMW(float *input, float *ref, float xp, float yp, float *ak, int order, float bg, int width, int height);
    
    
    
    void apply_planar_homography(float *input, int width, int height, laMatrix &H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight);
    void apply_planar_homography_zoom(float *input, int width, int height, laMatrix &H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight, float fZoom);
   
    
    void spline_interpolation_zoom(
                                   float *input,     // image to be warped
                                   const int    nx,        // image width
                                   const int    ny,        // image height
                                   const float  fFactor,   // zoom factor
                                   const int    order,     // order of interpolation
                                   const float bvalue,     // value outside the region
                                   float       *output     // image warped with bicubic interpolation
    );

    
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	//! Add Noise
	//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void compute_noise_sd(float **image, float *sd_noise, int width, int height, int num_channels, float SNR_db);
    void fpAddNoiseGaussian(float *u, float *v, float std, long int randinit, int size);
	void fpAddNoiseGaussianAfine(float *u,float *v, float a,float b,long int randinit, int size);
	
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	//! Gradient Computation
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void fiComputeImageGradient(float * fpI,float *fpXgrad, float *fpYgrad, float *fpGrad, float * fpOri, int iWidth, int iHeight, char cType = 'f');
    
	void fiComputeImageGradient(float * fpI, float *fpGrad, float * fpOri, int iWidth, int iHeight, char cType = 'f');
    
	void fiComputeImageGradient(float * fpI, float *fpGrad, int iWidth, int iHeight, char cType = 'f');
	
    
    
    
    
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Histogram
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
#define HSP_NB_CASES 500
    
    void fpHisto(float* input, laVector &histo, float *iminim, float *imaxim, int *n, float *s, int size, char flag);
    void fk_apply_histo(float *Histo,float *HistoInverse,int iRangeMax, float *src,float *srt, int width, int height);
    void fk_fill_histo(float *Histo,float *HistoInverse, int iRangeMax, float *in, int width, int height);
    void fk_histogram_specification(float *in1, float *in2, float *out, int width1, int height1, int width2, int height2);
    void fk_histogram_midway(float *in1, float *in2, float *out1, float *out2, int width1, int height1, int width2, int height2);
    void fk_histogram_midway_sequence(float **in, float **out, int nImages, int width, int height);
    void fk_apply_histo_mask(float *Histo,float *HistoInverse,int iRangeMax, float *src, float *mask, float *srt, int width, int height);
    void fk_fill_histo_mask(float *Histo,float *HistoInverse, int iRangeMax, float *in, float *mask, int width, int height);
    void fk_histogram_specification_mask(float *in1, float *mask1, float *in2, float *mask2, float *out, int width1, int height1, int width2, int height2);
    
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //! Geometrical Transformations
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void compute_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, laMatrix &H);
    int compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, laMatrix &H, int *accorded);
    int compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, laMatrix &H);
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //! Numerical Analysis Classes
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
  
    
    
    
    class laVector {
    protected:
        int d_n;	// size of array. upper index is nn-1
        float *d_v;
    public:
        
        laVector();
        explicit laVector(int n);								// Zero-based array
        laVector(const float &a, int n);						// initialize to constant value
        laVector(const float *a, int n);						// Initialize to array
        
        laVector(const laVector &rhs);							// Copy constructor
        laVector(const char * filename);
        
        laVector & operator=(const laVector &rhs);				// assignment
        laVector & operator=(const float &a);					// assign a to every element
        
		float & operator[](const int i);				// i'th element
		const float & operator[](const int i) const;
        
		
		
		float * v();
		
		int size() const;
        ~laVector();
        
        void erase();
        void create(int n);
        laVector copyBlock(int i0, int length);
        void sort(int decreasing_order);
        
        friend class laMatrix;
        
    };
    
    
    
    
	
    class laMatrix
    {
        
    protected:
        int d_n;			// nrows
        int d_m;			// ncols
        float **d_v;		// pointer
        
    public:
        
        
        //! Construction / Destruction
        
		laMatrix();
        laMatrix(int n, int m);
        laMatrix(const float &a, int n, int m);
        laMatrix(const float *a, int n, int m);
        laMatrix(const laMatrix &rhs);
        
        laMatrix & operator=(const laMatrix &rhs);
        laMatrix & operator=(const float &a);
        
        ~laMatrix();
        
        
		
        //! Basic operators
        
        float * operator[](const int i);	//subscripting: pointer to row i
        const float * operator[](const int i) const;
        
        float ** v();
        
        int nrows() const;
        int ncols() const;
        
        
        void create(int n, int m);
        
        
        
        //! Non member Arithmetic Operations
		
		friend laMatrix operator*  (float a, const laMatrix& rhs);                               // scalar matrix product
		friend laMatrix operator/  (const laMatrix& lhs, float a);                               // matrix scalar division
		friend laMatrix operator+  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix sum
		friend laMatrix operator-  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix subtraction
		friend laMatrix operator*  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix product
		
		friend laVector operator*  (const laMatrix& lhs, const laVector & rhs);                  // matrix vector product
		
        
		
		//! Other
		laMatrix transposed();
        
        laMatrix copyBlock(int i0, int j0, int rowb, int colb);
        
        friend class laVector;
        
    };
    
    
    
	void luinv(laMatrix &a, laMatrix &inv);
    void lusolve(laMatrix &a, laVector &x, laVector &b);
    void invCC(laMatrix &inv);          //! Inversion of definite positive matrices
    void compute_svd(laMatrix &A, laMatrix &m_U, laMatrix &m_V, laVector &m_W);
    void compute_svd_double(laMatrix &Ain, laMatrix &m_Uin, laMatrix &m_Vin, laVector &m_Win);
    void compute_pca_svd(laMatrix &X, laVector &S, laMatrix &V, laMatrix &U);
    
    
    void  l2_baricenter(float **X,float *baricenter,int n, int p);
    void  l1_baricenter(float **X,float *baricenter,int n, int p, int niter);
    
    
    void kmeans(float **data, float **means, int *labels, int nIter, int kCenters, int nVectors, int pCoordinates);
	void linear_fitting(float * tpX,float* tpY, float &da, float &db, int ilength);
    
    
    
    
    void  center_data_columns(float **X,float *baricenter,int n, int p);
    void  estimate_noise_pca_variances(int d, int n, float fSigma, float fDelta, int *osize, float *curve);
    
    
    
    //!
    //! Deprecated Numerical Analysis
    //!
    
	float ** allocate_float_matrix(int nrows, int ncols);
	void desallocate_float_matrix(float **matrix, int nrows, int ncols);
	
	
	/*- Householder reduction of a real symmetric matrix A[0..n-1][0..n_1]. On output, A is ---*/
	/*- replaced by the ortogonal matrix Q effecting the transformation. d[0..n-1] returns ----*/
	/*- the diagonal elements of the diagonal matrix, and e[0..n-1] the off-diagonal elements -*/
	/*- with e[0]=0. 									   */
	void symmetric_2_tridiag_HH(float **a,float *d,float *e,int n);
	
    
    /* QL with implicit shifts, to determine eigenvalues and eigenvectors of  a real, symmetric, tridiagonal matrix */
	/* d[0..n-1] contains the diagonal elements of the matrix and as output the returns the eigenvalues             */
	/* e[0..n-1] contains the sub-diagonal elements with e[0] arbitrary, on output e is destroyed			*/
	/* z[0..n-1][0..n-1] contain the identity matrix in input or the output of symmetric_2_tridiag_HH if previously */
	/* applied. 													*/
	
	int eigenvalues_QLI(float * d, float *e,float **z, int n);
	
    
    
	//void  pca_center_data(float **X,float *baricenter,int n, int p);
	//void order_decreasing(float *values, int *indexos, int size);
    
	float **  covariance_matrix(float **x,int n, int p);
    int  compute_pca_from_covariance(float **Mat_cov,float **Pcs, float *sVar, int p);
    
    
    
    

}



//
//! Parser
//

//! structure for parameters and options which are
//! optional in the code or they already have a default value
typedef struct optstruct
{
    char *gp;           //! string of two letters  "a:"  as necessary for using getopt afterwards
                        //! the ":" indicates that the activation of the option requires a value for the parameter
                        //! and "a" that this option is activated by "-a " in the command
    
    int flag;           //! flag indicating that the option has been activated
    
    char *defvalue;     //! default value for this parameter if not modified by console
    
    char *value;        //! value of the associated parameter to current option
    
    char *comment;      //! comment that appears by console
    
} OptStruct;



//! structure for necessary parameters of the method
typedef struct parstruct
{
    char * name;
    char * value;       //! value of the parameter
    char * comment;     //! comment that appears by console
    
} ParStruct;



int parsecmdline(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct*> & par);





namespace libUSTGDOUBLE
{

struct kernel_item {
        double v;
        int i, j;
    };


struct sorted_kernel {
        int nitems;
        int w, h;
        int ncut;
        struct kernel_item *items;
    };

	#ifndef PI
	#define PI 3.14159265358979323846264338327
	#endif

	#ifndef M_PI
	#define M_PI 3.14159265358979323846264338327
	#endif

	#define dTiny 1e-10
	#define fTiny 0.0000000001f
	#define fLarge 10000000000.0f
	#define dLarge 1e+10

	#define MAX(i,j) ( (i)<(j) ? (j):(i) )
	#define MIN(i,j) ( (i)<(j) ? (i):(j) )

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Double Value operations
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void fpClear(double *fpI,double fValue, int iLength);
	void fpCopy(double *fpI,double *fpO, int iLength);
	void fiPatchStatistics(double *fpIn, double *fpMinV, double *fpMaxV, double *fpMeanV, double *fpVarV, double *fpMedianV, double fRadius, int iWidth, int iHeight);
	void fiPatchMean(double *fpIn, double *fpMeanV, double fRadius, int iWidth, int iHeight);
	void fiPatchMax(double *fpIn, double *fpMaxV, double fRadius, int iWidth, int iHeight);
	void fiPatchMin(double *fpIn, double *fpMinV, double fRadius, int iWidth, int iHeight);
	void fiPatchVar(double *fpIn, double *fpVarV, double fRadius, int iWidth, int iHeight);
	void fiPatchMedian(double *fpIn, double *fpMedianV, double fRadius, int iWidth, int iHeight);
	void fiComputeImageGradient(double * fpI, double *fpGrad, double * fpOri, int iWidth, int iHeight, char cType);
	void fiComputeImageGradient(double * fpI, double *fpGrad, int iWidth, int iHeight, char cType);
	void fiComputeImageGradient(double * fpI,double *fpXgrad, double *fpYgrad, double *fpGrad, double * fpOri, int iWidth, int iHeight, char cType);
	void fiImageSampleAglomeration(double *igray,double *ogray, int factor, int width, int height);
	double*  fiFloatGaussKernel(double std, int & size);
	double*  fiDoubleGaussKernel(double std, int & size);
	void fiDoubleHorizontalConvolution(double *u, double *v, int width, int height, double *kernel, int ksize, int boundary);
	void fiDoubleVerticalConvolution(double *u, double *v, int width, int height, double *kernel,int ksize, int boundary);
	void fiConvol(double *u,double *v,int width,int height,double *kernel,int kwidth,int kheight, int boundary);
	int neumann_bc(int x, int nx, bool *out);
	int symmetric_bc(int x, int nx, bool *out);
	void fiConvol_skernel(double *u,double *v,int width,int height, struct sorted_kernel *skernel, int boundary);
	void fpAddNoiseGaussian(double *u, double *v, double std, long int randinit, int size);
	void fpAddNoiseGaussianAfine(double *u,double *v, double a,double b,long int randinit, int size);
	void nn_interpolation_zoom(
				 const double *input,     // image to be warped
				 const int    nx,        // image width
				 const int    ny,        // image height
				 const double  fFactor,     // zoom factor
				 double       *output    // image warped with bicubic interpolation
		 );

	void bicubic_interpolation_zoom(
			const double *input,     // image to be warped
			const int    nx,        // image width
			const int    ny,        // image height
			const double  fFactor,     // zoom factor
			const int    binterpolation,  // interpolate boundary if activated, else return bvalue
			const double bvalue,           // value outside the region
			double       *output    // image warped with bicubic interpolation
	);

	void bicubic_interpolation_zoom(
			const double *input,     // image to be warped
			const int    nx,        // image width
			const int    ny,        // image height
			const double  fFactor,     // zoom factor
			const int    binterpolation,  // interpolate boundary if activated, else return bvalue
			const double bvalue,           // value outside the region
			double       *output    // image warped with bicubic interpolation
	);

	void spline_interpolation_zoom(
			double *input,     // image to be warped
			const int    nx,        // image width
			const int    ny,        // image height
			const double  fFactor,   // zoom factor
			const int    order,     // order of interpolation
			const double bvalue,     // value outside the region
			double       *output     // image warped with bicubic interpolation
	);

	void finvsplineMW(double *in,int order,double *out, int width, int height);
	void init_splinenMW(double *a,int n);
	double evaluate_splineMW(double *input, double *ref, double xp, double yp, double *ak, int order, double bg, int width, int height);
	void invspline1DMW(double *c,int size,double *z,int npoles);
	void keysMW(double *c,double t,double a);
	void spline3MW(double *c,double t);
	void init_splinenMW(double *a,int n);
	 double vMW(double *in,int x,int y,double bg, int width, int height);
	 void splinenMW(double *c,double t,double *a,int n);
	 double initcausalMW(double *c,int n,double z);
	 double initanticausalMW(double *c,int n,double z);
	 double ipowMW(double x,int n);
	 void fiFFTZoom(double *in, double *out, double zoomx, double zoomy, int width, int height);
	 void fiFFTShearPascal(double dAmountOfShear, double *pFloatImageInput,double *pFloatImageOutput, int iAxis, double fDelta, int PleaseDoNotCenter, int width, int height);






	double fpMax(double *u,int *pos, int size);
	double fpMin(double *u,int *pos,int size);
	double fpVar(double *u,int size);
	double fpMean(double *u,int size);
	double fpMedian(double *u,int size);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Float pointer ordering
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void fpQuickSort(double *fpI, int iLength, int inverse = 0 );
	void fpQuickSort(double *fpI, double *fpO, int iLength, int inverse = 0);



	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//! Image Conversion
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	#define COEFF_YR 0.299
	#define COEFF_YG 0.587
	#define COEFF_YB 0.114

	void fiRgb2YuvDOUBLE(double *r,double *g,double *b,double *y,double *u,double *v,int width,int height);
	void fiRgb2YuvODOUBLE(double *r,double *g,double *b,double *y,double *u,double *v,int width,int height);
	void fiYuvO2RgbDOUBLE(double *r,double *g,double *b,double *y,double *u,double *v, int width,int height);
    void fiYuv2RgbDOUBLE(double *r,double *g,double *b,double *y,double *u,double *v, int width,int height);




}








#endif














