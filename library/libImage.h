#ifndef __libImage__
#define __libImage__



extern "C" {
#include "libImageFormats.h"
}


#define USE_FFTW
#ifdef USE_FFTW
#include <fftw3.h>
#endif

#include "libImageFormatPM.h"
#include "libBasic.h"


#define iipHorizontal 0
#define iipVertical 1


namespace libUSTG
{


//
//! Class definitions
//

class cflimage;
class flimage;


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Image Classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


class cflimage {


protected:



public:

	int d_c, d_w, d_h, d_wh, d_whc;		// channels, width, height, width*height, width*height*channels
	float *d_v;							// pointer


	char name[128];
	float visuMin, visuMax;				// useful only for visualization with mxview



	//
	//! Construction / Memory
	//

	cflimage();
	cflimage(int w, int h, int c);
	cflimage(int w, int h, float *igray);
	cflimage(int w, int h, float *ired, float *igreen, float *iblue);

	cflimage(int w, int h, unsigned char *igray);
	cflimage(int w, int h, unsigned char *ired, unsigned char *igreen, unsigned char *iblue);

	cflimage(const flimage &red, const flimage &green, const flimage &blue);


	void create(int w, int h, int c);

	cflimage(const cflimage& im);

	void erase();
	virtual ~cflimage();



	//
	//! Operators
	//

	cflimage& operator= (const cflimage& im);



	//
	//! Load / Save
	//
	void load(const char *filename);
	void save(const char *filename);



	//
	//! Get Basic Data
	//

	int c() const {return d_c;}
	int w() const {return d_w;}
	int h() const {return d_h;}
	int wh() const {return d_wh;}
	int whc() const {return d_whc;}


	float * v() {return d_v;}
	float * v(int i) {/*assert(i>= 0 && i < d_c);*/ return &d_v[i * d_wh];}

	inline float operator[](int i) const { /*assert(i>=0 && i < d_whc);*/ return d_v[i]; }
	inline float& operator[](int i) { /*assert(i>=0 && i < d_whc);*/ return d_v[i]; }


	flimage getChannel(int i);
	void getChannel(int i, flimage *out);

	operator  flimage();



	int isSameSize(const  cflimage &inIm);

	char * getName(){return name;}
	void  setName(const char *iname){strcpy(name, iname);}


	float getVisuMin(){return visuMin;}
	float getVisuMax(){return visuMax;}

	void setVisuMin(float a){visuMin = a;}
	void setVisuMax(float a){visuMax = a;}


	//
	//! Math
	//

	cflimage& operator= (float a);
	void operator-=(float a);
	void operator+=(float a);
	void operator*=(float a);

	float  min();
	float  max();

	float  min_channel (int i);
	float  max_channel (int i);

	void  max (float M);
	void  min (float m);

	void normalizeL1();

	void rint();
	void abs();

	void thre(float m, float M);



	//
	//! Color Conversion
	//

	flimage getGray();
	flimage getGray(float wr, float wg, float wb);
	flimage getGray(float wr, float wg, float wb, float wnir);

	void getGray(flimage * out);
	void getGray(float wr, float wg, float wb, flimage *out);
	void getGray(float wr, float wg, float wb, float wnir, flimage *out);


	void Rgb2Yuv(int iflagOrto);
	void Yuv2Rgb(int iflagOrto);


	cflimage binarize(float value, int inverse);
	void binarize(float value, int inverse, cflimage *out);


	//
	//! Block operations
	//

	cflimage copy(int ipx, int ipy, int iw, int ih);
	void paste(const cflimage &im, int x, int y);

	cflimage padding(int w, int h, float fValue);
	cflimage append(const cflimage &imIn, int extension);



	//
	//! Block operations
	//

	cflimage  patchMean(float fRadius);
	cflimage  patchVar(float fRadius);
	cflimage  patchMin(float fRadius);
	cflimage  patchMax(float fRadius);
	cflimage  patchMedian(float fRadius);


	cflimage  patchMean(flimage &kernel);
	cflimage  patchVar(flimage &kernel);



	//
	//! Gradient
	//
	cflimage gradient(char cType);
	cflimage xgradient(char cType);
	cflimage ygradient(char cType);
	cflimage gradient(cflimage &orientation, char cType);
	cflimage gradient(cflimage &xgrad, cflimage &ygrad, cflimage &orientation, char cType);




	//
	//! Patch distances
	//

	friend float distanceL2(const cflimage &input1,const  cflimage &input2, int ipx, int ipy, int iqx, int iqy);
	friend float distancePatchL2_NN(const  cflimage &input1,const  cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h);                         //! Non normalized distance
	friend float distancePatchL2_NN_Thresh(const  cflimage &input1,const  cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h, float fThresh);   //! Non normalized distance with threshold


	friend float distancePatchL2(const  cflimage &input1,const  cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h);


	friend float distancePatchL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h, cflimage &mean1, cflimage &mean2);


	friend float distancePatchWL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel);
	friend float distancePatchWL2_NN(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel);
	friend float distancePatchWL2_NN_Thresh(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel, float fThresh);


	friend float distancePatchWL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2);

	friend float distancePatchWL2M_NN_Rob(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2,  float fMinDist);

	friend float distancePatchWL2M_NN(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2);

	friend float distancePatchWL2M_NN_Thresh(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2, float fThresh);

	friend float distancePatchWL2_NN_Thresh_Rob(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel, float fThresh, float fMinDist);


	friend float distancePatchWL2M_NN_Thresh_Rob(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2, float fThresh, float fMinDist);


	//
	//! Subsampling & zooming & convolution
	//
	void convolveGauss(float fSigma, int boundary, cflimage *out);
	void convolve(const flimage &kernel, int boundary, cflimage *out);
	void convolve_skernel(struct sorted_kernel *skernel, int boundary, cflimage *out);


	cflimage subSample(int fFactor);
	void subSample(int fFactor, cflimage *out);


	cflimage subSampleAglomeration(int iFactor);
	void subSampleAglomeration(int iFactor, cflimage *out);

	cflimage subSampleConv(int fFactor, int boundary, float fSigma);
	void subSampleConv(int fFactor, int boundary, float fSigma, cflimage *out);


	void upSampleNN(const float fFactor, cflimage *out);

	void upSampleCubicSplines(const float fFactor, const int bintFlag, const float bValue, cflimage *out);
	cflimage upSampleCubicSplines(const float fFactor, const int bintFlag, const float bValue);

	void upSampleGenericSplines(const float fFactor, const int order, const float bValue, cflimage *out);


	cflimage mirror(int Orientation);


	cflimage UpSampleFFT(float fFactorx);
	cflimage UpSampleFFT(float fFactorx, float fFactory);


	void fftRot(float angle, float xtrans , float ytrans, int flagNoCenter, int flagSymmetric);
	void fftTrans(float xtrans,float ytrans, int flagNoCenter, int flagSymmetric);
	void fftShear(float dAmountOfShear, int iOrientation, float fDelta, int flagNoCenter, int flagSymmetric);





	//
	//! Value operations
	//
	void addGaussianNoise(float std);
	void addGaussianNoiseSD(float std, float sdstd);



};



class flimage : public cflimage
{

public:


	flimage();
	flimage(int w, int h);
	flimage(int w, int h, float *ptr);
	flimage(int w, int h, unsigned char *ptr);

	flimage(const flimage& im);

	void load(const char* filename);


	void create(int w, int h);


	using cflimage::operator=;
	flimage& operator=(const flimage& im);


	float operator()(int i, int j) const {return this->d_v[j * this->d_w + i];}
	float& operator()(int i, int j) { return this->d_v[j * this->d_w + i];}



	virtual ~flimage() {};

};









////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Movie Classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class cflmovie
{

public:


	// Constructor && Destructor
	cflmovie();
	cflmovie(const char * ifilename);							/// Reads a cflmovie file
	cflmovie(const char * ifilename,  int inframes);   			/// Opens a cflmovie file in writable format


	cflmovie(const cflmovie & imovie);

	~cflmovie();


	cflmovie& operator= (const cflmovie& im);


	// Load
	void load(const char * ifilename);											/// Reads a cflmovie file


	// Main variables
	int n() {return d_n;}



	// Getting image at position fpos
	cflimage getframe(int fpos);


	// Write image into film
	void write(cflimage &frame);


private:

	int d_n;

	char* filename;
	char* fileadress;				/// Used for reading/writting image names in format .mov

	bool  writable;				/// Is current cflmovie writable
	char ** strings;   			/// In case we read each time the image from the file

	int pos;        				/// Current position for writting
	std::ofstream outfile;     		/// Output mov file for writing

};





////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! 3D data
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class flData3D
{

public:


	///////////////// Constructor && Destructor
	flData3D();
	flData3D(flData3D &inData);
	flData3D(int inpoints, float *idx, float *idy, float *idz, float *idr, float *idg, float *idb);
	flData3D(const char * filename);

	~flData3D();


	///////////////// Operations
	flData3D& operator= (const flData3D &in);

	void loadFile(const char * filename);
	void SaveFile(const char * dataname);


	/////////////// Allocate Desallocate
	void allocate_coordinates(int n);
	void allocate_color(int n);
	void desallocate_coordinates();
	void desallocate_color();


	////////////// Get pointers
	float *getX(){return dx;}
	float *getY(){return dy;}
	float *getZ(){return dz;}

	float *getR(){return dr;}
	float *getG(){return dg;}
	float *getB(){return db;}


	///////////// Get number of points
	int getSize(){return npoints;}


	bool Ok(){return ok;}
	bool Color(){return color;}



	void normalize();
	//void normalize(float bx, float by, float bz, float amp);
	//void getNormParam(float *bx, float *by,float *bz, float *amp);



public:

	int npoints;
	float *dx, *dy, *dz;
	float *dr, *dg, *db;

	bool ok;
	bool color;

};








}


namespace libUSTGDOUBLE
{

//
//! Class definitions
//

class cflimageDOUBLE;
class flimageDOUBLE;


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Image Classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class cflimageDOUBLE {


protected:



public:

	int d_c, d_w, d_h, d_wh, d_whc;		// channels, width, height, width*height, width*height*channels
	double *d_v;							// pointer


	char name[128];
	double visuMin, visuMax;				// useful only for visualization with mxview



	//
	//! Construction / Memory
	//

	cflimageDOUBLE();
	cflimageDOUBLE(int w, int h, int c);
	cflimageDOUBLE(int w, int h, double *igray);
	cflimageDOUBLE(int w, int h, double *ired, double *igreen, double *iblue);

	cflimageDOUBLE(int w, int h, unsigned char *igray);
	cflimageDOUBLE(int w, int h, unsigned char *ired, unsigned char *igreen, unsigned char *iblue);

	cflimageDOUBLE(const flimageDOUBLE &red, const flimageDOUBLE &green, const flimageDOUBLE &blue);


	void create(int w, int h, int c);

	cflimageDOUBLE(const cflimageDOUBLE& im);

	void erase();
	virtual ~cflimageDOUBLE();



	//
	//! Operators
	//

	cflimageDOUBLE& operator= (const cflimageDOUBLE& im);



	//
	//! Load / Save
	//
	void load(const char *filename);
	void save(const char *filename);



	//
	//! Get Basic Data
	//

	int c() const {return d_c;}
	int w() const {return d_w;}
	int h() const {return d_h;}
	int wh() const {return d_wh;}
	int whc() const {return d_whc;}


	double * v() {return d_v;}
	double * v(int i) {/*assert(i>= 0 && i < d_c);*/ return &d_v[i * d_wh];}

	inline double operator[](int i) const { /*assert(i>=0 && i < d_whc);*/ return d_v[i]; }
	inline double& operator[](int i) { /*assert(i>=0 && i < d_whc);*/ return d_v[i]; }


	flimageDOUBLE getChannel(int i);
	void getChannel(int i, flimageDOUBLE *out);

	operator  flimageDOUBLE();



	int isSameSize(const  cflimageDOUBLE &inIm);

	char * getName(){return name;}
	void  setName(const char *iname){strcpy(name, iname);}


	double getVisuMin(){return visuMin;}
	double getVisuMax(){return visuMax;}

	void setVisuMin(double a){visuMin = a;}
	void setVisuMax(double a){visuMax = a;}


	//
	//! Math
	//

	cflimageDOUBLE& operator= (double a);
	void operator-=(double a);
	void operator+=(double a);
	void operator*=(double a);

	double  min();
	double  max();

	double  min_channel (int i);
	double  max_channel (int i);

	void  max (double M);
	void  min (double m);

	void normalizeL1();

	void rint();
	void abs();

	void thre(double m, double M);



	//
	//! Color Conversion
	//

	flimageDOUBLE getGray();
	flimageDOUBLE getGray(double wr, double wg, double wb);
	flimageDOUBLE getGray(double wr, double wg, double wb, double wnir);

	void getGray(flimageDOUBLE * out);
	void getGray(double wr, double wg, double wb, flimageDOUBLE *out);
	void getGray(double wr, double wg, double wb, double wnir, flimageDOUBLE *out);


	void Rgb2YuvDOUBLE(int iflagOrto);
	void Yuv2RgbDOUBLE(int iflagOrto);


	cflimageDOUBLE binarize(double value, int inverse);
	void binarize(double value, int inverse, cflimageDOUBLE *out);


	//
	//! Block operations
	//

	cflimageDOUBLE copy(int ipx, int ipy, int iw, int ih);
	void paste(const cflimageDOUBLE &im, int x, int y);

	cflimageDOUBLE padding(int w, int h, double dValue);
	cflimageDOUBLE append(const cflimageDOUBLE &imIn, int extension);



	//
	//! Block operations
	//

	cflimageDOUBLE  patchMean(double fRadius);
	cflimageDOUBLE  patchVar(double fRadius);
	cflimageDOUBLE  patchMin(double fRadius);
	cflimageDOUBLE  patchMax(double fRadius);
	cflimageDOUBLE  patchMedian(double fRadius);


	cflimageDOUBLE  patchMean(flimageDOUBLE &kernel);
	cflimageDOUBLE  patchVar(flimageDOUBLE &kernel);



	//
	//! Gradient
	//
	cflimageDOUBLE gradient(char cType);
	cflimageDOUBLE xgradient(char cType);
	cflimageDOUBLE ygradient(char cType);
	cflimageDOUBLE gradient(cflimageDOUBLE &orientation, char cType);
	cflimageDOUBLE gradient(cflimageDOUBLE &xgrad, cflimageDOUBLE &ygrad, cflimageDOUBLE &orientation, char cType);




	//
	//! Patch distances
	//

	friend double distanceL2(const cflimageDOUBLE &input1,const  cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy);
	friend double distancePatchL2_NN(const  cflimageDOUBLE &input1,const  cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h);                         //! Non normalized distance
	friend double distancePatchL2_NN_Thresh(const  cflimageDOUBLE &input1,const  cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h, double fThresh);   //! Non normalized distance with threshold


	friend double distancePatchL2(const  cflimageDOUBLE &input1,const  cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h);


	friend double distancePatchL2M(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h, cflimageDOUBLE &mean1, cflimageDOUBLE &mean2);


	friend double distancePatchWL2(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel);
	friend double distancePatchWL2_NN(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel);
	friend double distancePatchWL2_NN_Thresh(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel, double fThresh);


	friend double distancePatchWL2M(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel,  cflimageDOUBLE &mean1, cflimageDOUBLE &mean2);

	friend double distancePatchWL2M_NN_Rob(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel,  cflimageDOUBLE &mean1, cflimageDOUBLE &mean2,  double fMinDist);

	friend double distancePatchWL2M_NN(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel,  cflimageDOUBLE &mean1, cflimageDOUBLE &mean2);

	friend double distancePatchWL2M_NN_Thresh(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel,  cflimageDOUBLE &mean1, cflimageDOUBLE &mean2, double fThresh);

	friend double distancePatchWL2_NN_Thresh_Rob(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel, double fThresh, double fMinDist);


	friend double distancePatchWL2M_NN_Thresh_Rob(cflimageDOUBLE &input1, cflimageDOUBLE &input2, int ipx, int ipy, int iqx, int iqy, flimageDOUBLE &kernel,  cflimageDOUBLE &mean1, cflimageDOUBLE &mean2, double fThresh, double fMinDist);


	//
	//! Subsampling & zooming & convolution
	//
	void convolveGauss(double fSigma, int boundary, cflimageDOUBLE *out);
	void convolve(const flimageDOUBLE &kernel, int boundary, cflimageDOUBLE *out);
	void convolve_skernel(struct sorted_kernel *skernel, int boundary, cflimageDOUBLE *out);


	cflimageDOUBLE subSample(int fFactor);
	void subSample(int fFactor, cflimageDOUBLE *out);


	cflimageDOUBLE subSampleAglomeration(int iFactor);
	void subSampleAglomeration(int iFactor, cflimageDOUBLE *out);

	cflimageDOUBLE subSampleConv(int fFactor, int boundary, double fSigma);
	void subSampleConv(int fFactor, int boundary, double fSigma, cflimageDOUBLE *out);


	void upSampleNN(const double fFactor, cflimageDOUBLE *out);

	void upSampleCubicSplines(const double fFactor, const int bintFlag, const double bValue, cflimageDOUBLE *out);
	cflimageDOUBLE upSampleCubicSplines(const double fFactor, const int bintFlag, const double bValue);

	void upSampleGenericSplines(const double fFactor, const int order, const double bValue, cflimageDOUBLE *out);


	cflimageDOUBLE mirror(int Orientation);


	cflimageDOUBLE UpSampleFFT(double dFactorx);
	cflimageDOUBLE UpSampleFFT(double dFactorx, double dFactory);


	void fftRot(double angle, double xtrans , double ytrans, int flagNoCenter, int flagSymmetric);
	void fftTrans(double xtrans,double ytrans, int flagNoCenter, int flagSymmetric);
	void fftShear(double dAmountOfShear, int iOrientation, double fDelta, int flagNoCenter, int flagSymmetric);





	//
	//! Value operations
	//
	void addGaussianNoise(double std);
	void addGaussianNoiseSD(double std, double sdstd);



};


class flimageDOUBLE : public cflimageDOUBLE{

public:


	flimageDOUBLE();
	flimageDOUBLE(int w, int h);
	flimageDOUBLE(int w, int h, double *ptr);
	flimageDOUBLE(int w, int h, unsigned char *ptr);

	flimageDOUBLE(const flimageDOUBLE& im);

	void load(const char* filename);


	void create(int w, int h);


	using cflimageDOUBLE::operator=;
	flimageDOUBLE& operator=(const flimageDOUBLE& im);


	double operator()(int i, int j) const {return this->d_v[j * this->d_w + i];}
	double& operator()(int i, int j) { return this->d_v[j * this->d_w + i];}


	virtual ~flimageDOUBLE() {};

};

}

#endif
