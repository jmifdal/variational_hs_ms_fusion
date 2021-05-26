#ifndef _LIBFILTERINGVARIATIONAL_H
#define _LIBFILTERINGVARIATIONAL_H

#include "libBasic.h"


namespace libUSTGFiltVar
{
    
    /***** ROUTINES RELATED TO PRIMAL-DUAL ALGORITHM *****/
    
    // Filtering an image by the variational model composed of TV + NL regularizations and weighted L2 fidelity term.
    // Primal-dual hybrid gradient method is used.
    void filtering_TV_NL_L2w(float *u, float *f, float *mask, std::vector< std::vector<float> > wxy, std::vector< std::vector<int> > posxy,
                             std::vector< std::vector<float> > wyx, std::vector< std::vector<int> > posyx, std::vector< std::vector<int> > posw,
                             float lmbTV, float lmbNL, float tau, float sigma, float tol, int maxIter, int width, int height);
    
    // Compute proximity operator of L2 data-fitting term
    void proxGL2(float *u_upd, float *u, float *f, float *mask, float lmbFit, float tau,
               float *div, float *nldiv, int dim);
    
    
    // Compute proximity operator of L1 data-fitting term
    void proxGL1(float *u_upd, float *u, float *f, float *mask, float lmbFit, float tau,
                 float *div, float *nldiv, int dim);
    
    // Compute proximity operator of TV regularization
    void proxTV(float *px_upd, float *py_upd, float *px, float *py, float *gradx_upd, float *grady_upd,
                float *gradx, float *grady, float lmbTV, float sigma, int dim);
    
    // Compute proximity operator of NL regularization
    void proxNL(std::vector< std::vector <float> > &q_upd, std::vector< std::vector <float> > &q,
                std::vector< std::vector <float> > &nlgrad_upd, std::vector< std::vector <float> > &nlgrad,
                float lmbNL, float sigma, int dim);
    
    // Compute local gradient at each pixel via forward differences
    void forward_gradient(float *f, float *fx, float *fy, int nx, int ny);
    
    // Compute local divergence at each pixel as minus the adjoint operator of the local gradient, div = -nabla^T
    void divergence(float *v1, float *v2, float *div, int nx, int ny);
    
    // Compute nonlocal gradient operator at each pixel
    void nl_gradient(float *data, std::vector< std::vector<float> > &nlgrad, std::vector< std::vector<float> > &weights,
                     std::vector< std::vector<int> > &positions, int dim);
    
    // Compute nonlocal divergence at each pixel as minus the adjoint operator of the nonlocal gradient, div_w = -nabla_w^T
    void nl_divergence(std::vector< std::vector<float> > &data, float *nldiv, std::vector< std::vector<float> > &wxy,
                       std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx,
                       std::vector< std::vector<int> > &posw, int dim);
    
    // Compute error between two consecutive iterations
    float relative_error(float *u, float *u_upd, int dim);
    
    
    
    /***** ROUTINES RELATED TO NONLOCAL REGULARIZATION *****/

    // Compute weights for NL regularization on a prescribed image using all neighbouring pixels
    void nlweights_reg(float **image, std::vector< std::vector<float> > &wxy, std::vector< std::vector<int> > &posxy,
                       std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx,
                       std::vector< std::vector<int> > &posw, float hSim, float hClose, int numNeighbours,
                       int reswind, int compwind, int flagNormalize, int num_channels, int width, int height);
    
    // Compute weights for NL regularization on a prescribed image using all neighbouring pixels
    void nlweights_reg_mask(float **image, float *mask, std::vector< std::vector<float> > &wxy,
                            std::vector< std::vector<int> > &posxy, std::vector< std::vector<float> > &wyx,
                            std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw,
                            float hSim, float hClose, int numNeighbours, int reswind, int compwind, int flagNormalize,
                            int num_channels, int width, int height);
    
} // libUSTGFiltVar

#endif
