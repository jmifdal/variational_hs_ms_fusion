#ifndef _LIBHYPERSPECTRAL_H_
#define _LIBHYPERSPECTRAL_H_

//#include "libBasic.h"
#include "./libBasic.h"


namespace libUSTGHYPER
{

    // NL regularization with different weights for each hyperspectral channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L1 data constraint.
    void hyperfusion_NL_L1_wh(float **u_upd, float **fH, float **fM, float **S, float **St, float **PH, float **PHint, float **fHint,
                              std::vector< std::vector< std::vector<float> > > wxy, std::vector< std::vector< std::vector<int> > > posxy,
                              std::vector< std::vector< std::vector<float> > > wyx, std::vector< std::vector< std::vector<int> > > posyx,
                              std::vector< std::vector< std::vector<int> > > posw, float lmbH, float lmbM, float mu, float tau, float sigma,
                              float tol, int maxIter, int sampling, float stdBlur, int hs_channels, int ms_channels, int width, int height);
	//DEBUGG
	/*void hyperfusion_NL_L1_wh(float **u_upd, float **fH, float **fM, float **S, float **St, float **PH, float **PHint, float **fHint,
                              std::vector< std::vector< std::vector<float> > > wxy, std::vector< std::vector< std::vector<int> > > posxy,
                              std::vector< std::vector< std::vector<float> > > wyx, std::vector< std::vector< std::vector<int> > > posyx,
                              std::vector< std::vector< std::vector<int> > > posw, float lmbH, float lmbM, float mu, float tau, float sigma,
                              float tol, int maxIter, int sampling, float stdBlur, int hs_channels, int ms_channels, int width, int height,
                              std::vector< std::vector< std::vector<float> > > wxy_L2, std::vector< std::vector< std::vector<int> > > posxy_L2,
                              std::vector< std::vector< std::vector<float> > > wyx_L2, std::vector< std::vector< std::vector<int> > > posyx_L2,
                              std::vector< std::vector< std::vector<int> > > posw_L2);*/

    
    // NL regularization with the same weight for all hyperspectral channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L1 data constraint.
    void hyperfusion_NL_L1(float **u_upd, float **fH, float **fM, float **S, float **St, float **PH, float **PHint, float **fHint,
                           std::vector< std::vector<float> > wxy, std::vector< std::vector<int> > posxy, std::vector< std::vector<float> > wyx,
                           std::vector< std::vector<int> > posyx, std::vector< std::vector<int> > posw, float lmbH, float lmbM, float mu,
                           float tau, float sigma, float tol, int maxIter, int sampling, float stdBlur, int hs_channels, int ms_channels,
                           int width, int height);
    
    
    // NL regularization with different weights for each hyperspectral channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L2 data constraint.
    void hyperfusion_NL_wh(float **u_upd, float **fH, float **fM, float **S, float **St, float **PH, float **PHint, float **fHint,
                           std::vector< std::vector< std::vector<float> > > whxy, std::vector< std::vector< std::vector<int> > > poshxy,
                           std::vector< std::vector< std::vector<float> > > whyx, std::vector< std::vector< std::vector<int> > > poshyx,
                           std::vector< std::vector< std::vector<int> > > poshw, float lmbH, float lmbM, float mu, float tau, float sigma,
                           float tol, int maxIter, int sampling, float stdBlur, int hs_channels, int ms_channels, int width, int height);
    
    // NL regularization with the same weight for all hyperspectral channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L2 data constraint.
    void hyperfusion_NL(float **u_upd, float **fH, float **fM, float **S, float **St, float **PH, float **PHint, float **fHint,
                        std::vector< std::vector<float> > wxy, std::vector< std::vector<int> > posxy, std::vector< std::vector<float> > wyx,
                        std::vector< std::vector<int> > posyx, std::vector< std::vector<int> > posw, float lmbH, float lmbM, float mu,
                        float tau, float sigma, float tol, int maxIter, int sampling, float stdBlur, int hs_channels, int ms_channels,
                        int width, int height);
    
    // Compute proximity operator of dualized squared L2 data term
    void proxDualSquaredL2(float **p_upd, float **p_arg, float **f, float lmb, float sigma, int num_channels, int dim);
    
    // Compute proximity operator of NL regularization.
    void proxNL(std::vector< std::vector< std::vector <float> > > &p_upd, std::vector< std::vector< std::vector <float> > > &p_arg,
                int num_channels, int dim);
    
    // For each hyperspectral channel, compute the nonlocal gradient.
    // The weights are different for each hyperspectral channel.
    void nl_gradient(float **data, std::vector< std::vector< std::vector<float> > > &nlgrad,
                     std::vector< std::vector< std::vector<float> > > &wxy, std::vector< std::vector< std::vector<int> > > &posxy,
                     int num_channels, int dim);
    
    // For each hyperspectral channel, compute the nonlocal gradient.
    // The weights are the same for all hyperspectral channels.
    void nl_gradient(float **data, std::vector< std::vector< std::vector<float> > > &nlgrad, std::vector< std::vector<float> > &wxy,
                     std::vector< std::vector<int> > &posxy, int num_channels, int dim);
    
    
    // For each hyperspectral channel, compute the nonlocal divergence as div_w = -nabla_w^T.
    // The weights are different for each hyperspectral channel.
    void nl_divergence(std::vector< std::vector< std::vector<float> > > &data, float **nldiv,
                       std::vector< std::vector< std::vector<float> > > &wxy, std::vector< std::vector< std::vector<float> > > &wyx,
                       std::vector< std::vector< std::vector<int> > > &posyx, std::vector< std::vector< std::vector<int> > > &posw,
                       int num_channels, int dim);
    
    // For each hyperspectral channel, compute the nonlocal divergence as div_w = -nabla_w^T.
    // The weights are the same for all hyperspectral channels.
    void nl_divergence(std::vector< std::vector< std::vector<float> > > &data, float **nldiv, std::vector< std::vector<float> > &wxy,
                       std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx,
                       std::vector< std::vector<int> > &posw, int num_channels, int dim);
    
    // Compute Gaussian convolution in Fourier domain
    void FFT_gaussian_convol(float *convolved, float *data, float stdBlur, int width, int height);
    
    // Downsampling
    void downsampling(float *downsampled, float *data, int sampling_factor, int width, int height);
    
    // Upsampling by zero-padding
    void upsampling_zeropadding(float *upsampled, float *data, int sampling_factor, int width, int height);
    
    // Compute L2 error between two consecutive iterations
    float compute_error(float **u, float **u_upd, int num_channels, int dim);

} // libUSTGHYPER

#endif
