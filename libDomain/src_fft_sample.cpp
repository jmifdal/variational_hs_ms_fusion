#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;


int main(int argc, char **argv)
{
	
	
	vector <OptStruct *> options;
	OptStruct oA = {"a:", 0, "0.1",  NULL, "value at cutting frequency"};  options.push_back(&oA);	
	OptStruct oS = {"s:", 0, "2",  NULL, "sampling value"};  options.push_back(&oS);	
	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output image"}; parameters.push_back(&pout);
	
	
	
	if (!parsecmdline("src_fft_sample", "Samples with convolution of fixed value at cutting frequency", argc, argv, options, parameters))
		return 0;

	
	// parameters
	float sampling = atof(oS.value);
	float alpha = atof(oA.value);
	
	
	// input image
	libUSTG::cflimage input;
	input.load(pinput.value);
	
	if (input.w() != input.h())
	{
		printf("warning :: fft_sampling :: current version only deals with square images\n");

	}
	
	
	
	// kernel
	float sigma = (float) sqrt( - (4.0 * sampling * sampling) * log(alpha));
	sigma *= sigma; 
	
	float factor =  1.0  / (float) input.w();
	factor *= factor;

	
	libUSTG::flimage kernel0(input.w(), input.h());
	libUSTG::flimage kernel(input.w(), input.h());
	
	for(int j= - input.h()/2; j <= input.h() / 2 - 1; j++)
		for(int i= - input.w()/2; i <= input.w() / 2 -1 ; i++)
		{
			
			float dist = - factor * sigma * (double) (i*i + j*j);
			kernel0[input.w()/2 + i + (input.h()/2 + j) * kernel0.w()] = expf(dist);
			
		}
	
	
	for(int j= 0; j < input.h(); j++)
		for(int i= 0; i < input.w()  ; i++)
		{
			
			int ni = (i + input.w()/2) % input.w();
			int nj = (j + input.h()/2) % input.h();
			
			kernel[ni + nj*kernel.w()] = kernel0[i+ j * kernel0.w()];			
			
		}
	
	
	
	printf("debug :: fft_sampling :: kernel at cutting frequency %f\n",  kernel[input.w() / (2 * (int) sampling)]);

	
	
	
	// output
	libUSTG::cflimage output((int) floor( (double) input.w() / sampling),
				   (int) floor( (double) input.h() / sampling), 
				   input.c());
	
	
	for (int ii=0; ii < input.c(); ii++) 
	{

		libUSTG::flimage reel(input.w(), input.h());
		libUSTG::flimage imaginary(input.w(), input.h());
		libUSTG::flimage filtered(input.w(), input.h());
		libUSTG::flimage auxiliar(input.w(), input.h());

		// fft
		libUSTG::fft2d(input.v(ii), NULL,reel.v(),imaginary.v(),0,input.w(),input.h());
	
		
		// multiply by kernel
		for (int j=0; j < input.wh(); j++) 
		{
			reel[j] *= kernel[j];
			imaginary[j] *= kernel[j];
		}
		
		
		// ifft
		libUSTG::fft2d(reel.v(),imaginary.v(), filtered.v(),auxiliar.v(),1,input.w(),input.h());
	
		
		// sampling
		libUSTG::fiImageSample(filtered.v(), output.v(ii), sampling, input.w(), input.h());
		
	}
	
	
	output.save(pout.value); 
	
	
}



/*
 

## Sample images 

# Kernel takes value 0.1 at w/2/3
mxcnesresolution kernel 2052 9.1

fft2d -A reel -B imagin left10cm
fop -A kernel -t reel reel
fop -A kernel -t imagin imagin
fft2d -i -I imagin -A output reel

fsample output -ftype PM_F left30cm 3 

fft2d -A reel -B imagin right10cm
fop -A kernel -t reel reel
fop -A kernel -t imagin imagin
fft2d -i -I imagin -A output reel

fsample output -ftype PM_F right30cm 3

rm kernel reel imagin output




int main(int argc, char **argv)
{

		
	OptStruct *opt = NULL;

	ParStruct par[4] = {
				{"out", NULL, "output file"},
				{"width", NULL, "output width"},
				{"alpha", NULL, "value at cut frequency"},
				{"sigma", NULL, "sampling to be applied after"}
	};
	

	if (!parsecmdline("samplink_kernel","sampling_kernel", argc, argv, opt,0, par, 4))
		return 0;


	int width = atoi(par[1].value);
	int height = width;

	// sigma = sqrt(- (2 * sampling) ^ 2 sqrt(alpha))
	// being alpha the value of the kernel at coupure for a fixed sampling
	double alpha =  atof(par[2].value);
	double sampling = atof(par[3].value);

	double sigma = sqrt( - (4.0 * sampling * sampling) * log(alpha));
	sigma *= sigma; 
		
	double factor =  1.0  / (double) width;
	factor *= factor;


	float *igray = new float[width*height];
	for(int j= - height/2; j <= height / 2 - 1; j++)
		for(int i= - width/2; i <= width / 2 -1 ; i++)
		{


			double dist = - factor * sigma * (double) (i*i + j*j);
			igray[ (height/2 + j) * width + width/2 + i ] = exp(dist);
			

		}



	vflimage output(width, height , 1);
	float *ogray = output.getChannelPlane(0);

 
	for(int j= 0; j < height; j++)
		for(int i= 0; i < width  ; i++)
		{

			int ni = (i + width/2) % width;
			int nj = (j + height/2) % height;

			ogray[ nj*width + ni ] = igray[ j*width + i];			

		}

	output.save(par[0].value, "float"); 
	delete[] igray;

	printf("sigma: %f fe/s: %f\n", sqrt(sigma), ogray[width/(2 * (int) sampling)]);

	return 1;
}

*/