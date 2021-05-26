#ifndef _LIBHYPERSPECTRALDOUBLE_H_
#define _LIBHYPERSPECTRALDOUBLE_H_

//#include "libBasic.h"
#include "./libBasic.h"


namespace libUSTGHYPERDOUBLE
{

    // NL regularization with different weights for each hyperspectral channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L1 data constraint.
    void hyperfusion_NL_L1_wh(double **u_upd, double **fH, double **fM, double **S, double **St, double **PH, double **PHint, double **fHint,
                              std::vector< std::vector< std::vector<double> > > wxy, std::vector< std::vector< std::vector<int> > > posxy,
                              std::vector< std::vector< std::vector<double> > > wyx, std::vector< std::vector< std::vector<int> > > posyx,
                              std::vector< std::vector< std::vector<int> > > posw, double lmbH, double lmbM, double mu, double tau, double sigma,
                              double tol, int maxIter, int sampling, double stdBlur, int hs_channels, int ms_channels, int width, int height);
	//DEBUGG
	/*void hyperfusion_NL_L1_wh(double **u_upd, double **fH, double **fM, double **S, double **St, double **PH, double **PHint, double **fHint,
                              std::vector< std::vector< std::vector<double> > > wxy, std::vector< std::vector< std::vector<int> > > posxy,
                              std::vector< std::vector< std::vector<double> > > wyx, std::vector< std::vector< std::vector<int> > > posyx,
                              std::vector< std::vector< std::vector<int> > > posw, double lmbH, double lmbM, double mu, double tau, double sigma,
                              double tol, int maxIter, int sampling, double stdBlur, int hs_channels, int ms_channels, int width, int height,
                              std::vector< std::vector< std::vector<double> > > wxy_L2, std::vector< std::vector< std::vector<int> > > posxy_L2,
                              std::vector< std::vector< std::vector<double> > > wyx_L2, std::vector< std::vector< std::vector<int> > > posyx_L2,
                              std::vector< std::vector< std::vector<int> > > posw_L2);*/


    // NL regularization with the same weight for all hyperspectral channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L1 data constraint.
    void hyperfusion_NL_L1(double **u_upd, double **fH, double **fM, double **S, double **St, double **PH, double **PHint, double **fHint,
                           std::vector< std::vector<double> > wxy, std::vector< std::vector<int> > posxy, std::vector< std::vector<double> > wyx,
                           std::vector< std::vector<int> > posyx, std::vector< std::vector<int> > posw, double lmbH, double lmbM, double mu,
                           double tau, double sigma, double tol, int maxIter, int sampling, double stdBlur, int hs_channels, int ms_channels,
                           int width, int height);


    // NL regularization with different weights for each hyperspectral channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L2 data constraint.
    void hyperfusion_NL_wh(double **u_upd, double **fH, double **fM, double **S, double **St, double **PH, double **PHint, double **fHint,
                           std::vector< std::vector< std::vector<double> > > whxy, std::vector< std::vector< std::vector<int> > > poshxy,
                           std::vector< std::vector< std::vector<double> > > whyx, std::vector< std::vector< std::vector<int> > > poshyx,
                           std::vector< std::vector< std::vector<int> > > poshw, double lmbH, double lmbM, double mu, double tau, double sigma,
                           double tol, int maxIter, int sampling, double stdBlur, int hs_channels, int ms_channels, int width, int height);

    // NL regularization with the same weight for all hyperspectral channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    // L2 fidelity-term forcing closeness to hyperspectral data on the low-resolution domain, defined by low-pass filtering + subsampling.
    // L2 fidelity-term forcing closeness to multispectral data on the high-resolution domain, defined by matrix S.
    // L2 data constraint.
    void hyperfusion_NL(double **u_upd, double **fH, double **fM, double **S, double **St, double **PH, double **PHint, double **fHint,
                        std::vector< std::vector<double> > wxy, std::vector< std::vector<int> > posxy, std::vector< std::vector<double> > wyx,
                        std::vector< std::vector<int> > posyx, std::vector< std::vector<int> > posw, double lmbH, double lmbM, double mu,
                        double tau, double sigma, double tol, int maxIter, int sampling, double stdBlur, int hs_channels, int ms_channels,
                        int width, int height);

    // Compute proximity operator of dualized squared L2 data term
    void proxDualSquaredL2(double **p_upd, double **p_arg, double **f, double lmb, double sigma, int num_channels, int dim);

    // Compute proximity operator of NL regularization.
    void proxNL(std::vector< std::vector< std::vector <double> > > &p_upd, std::vector< std::vector< std::vector <double> > > &p_arg,
                int num_channels, int dim);

    // For each hyperspectral channel, compute the nonlocal gradient.
    // The weights are different for each hyperspectral channel.
    void nl_gradient(double **data, std::vector< std::vector< std::vector<double> > > &nlgrad,
                     std::vector< std::vector< std::vector<double> > > &wxy, std::vector< std::vector< std::vector<int> > > &posxy,
                     int num_channels, int dim);

    // For each hyperspectral channel, compute the nonlocal gradient.
    // The weights are the same for all hyperspectral channels.
    void nl_gradient(double **data, std::vector< std::vector< std::vector<double> > > &nlgrad, std::vector< std::vector<double> > &wxy,
                     std::vector< std::vector<int> > &posxy, int num_channels, int dim);


    // For each hyperspectral channel, compute the nonlocal divergence as div_w = -nabla_w^T.
    // The weights are different for each hyperspectral channel.
    void nl_divergence(std::vector< std::vector< std::vector<double> > > &data, double **nldiv,
                       std::vector< std::vector< std::vector<double> > > &wxy, std::vector< std::vector< std::vector<double> > > &wyx,
                       std::vector< std::vector< std::vector<int> > > &posyx, std::vector< std::vector< std::vector<int> > > &posw,
                       int num_channels, int dim);

    // For each hyperspectral channel, compute the nonlocal divergence as div_w = -nabla_w^T.
    // The weights are the same for all hyperspectral channels.
    void nl_divergence(std::vector< std::vector< std::vector<double> > > &data, double **nldiv, std::vector< std::vector<double> > &wxy,
                       std::vector< std::vector<double> > &wyx, std::vector< std::vector<int> > &posyx,
                       std::vector< std::vector<int> > &posw, int num_channels, int dim);

    // Compute Gaussian convolution in Fourier domain
    void FFT_gaussian_convol(double *convolved, double *data, double stdBlur, int width, int height);

    // Downsampling
    void downsampling(double *downsampled, double *data, int sampling_factor, int width, int height);

    // Upsampling by zero-padding
    void upsampling_zeropadding(double *upsampled, double *data, int sampling_factor, int width, int height);

    // Compute L2 error between two consecutive iterations
    double compute_error(double **u, double **u_upd, int num_channels, int dim);

} // libUSTGHYPERDOUBLE

#endif
