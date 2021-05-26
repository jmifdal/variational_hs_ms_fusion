#include "libFilteringVariational.h"

using namespace libUSTG;

namespace libUSTGFiltVar
{
    // Filtering an image by the variational model composed of TV + NL regularizations and weighted L2 fidelity term.
    // Primal-dual hybrid gradient method is used.
    void filtering_TV_NL_L2w(float *u, float *f, float *mask, std::vector< std::vector<float> > wxy, std::vector< std::vector<int> > posxy,
                             std::vector< std::vector<float> > wyx, std::vector< std::vector<int> > posyx, std::vector< std::vector<int> > posw,
                             float lmbTV, float lmbNL, float tau, float sigma, float tol, int maxIter, int width, int height)
    {
        
        // Image size
        int dim = width * height;
        
        // Primal variables
        float *u_upd = new float[dim];
        
        fpClear(u, 0.0f, dim);
        fpClear(u_upd, 0.0f, dim);
        
        // Dual variables related to TV regularization
        float *px = new float[dim];
        float *py = new float[dim];
        float *px_upd = new float[dim];
        float *py_upd = new float[dim];
        
        fpClear(px, 0.0f, dim);
        fpClear(py, 0.0f, dim);
        fpClear(px_upd, 0.0f, dim);
        fpClear(py_upd, 0.0f, dim);
        
        float *gradx = new float[dim];
        float *grady = new float[dim];
        float *gradx_upd = new float[dim];
        float *grady_upd = new float[dim];
        float *div = new float[dim];
        float *div_upd = new float[dim];
        
        fpClear(gradx, 0.0f, dim);
        fpClear(grady, 0.0f, dim);
        fpClear(gradx_upd, 0.0f, dim);
        fpClear(grady_upd, 0.0f, dim);
        fpClear(div, 0.0f, dim);
        fpClear(div_upd, 0.0f, dim);
        
        // Dual variables related to NLTV regularization
        std::vector< std::vector<float> > q;
        std::vector< std::vector<float> > q_upd;
        std::vector< std::vector<float> > nlgrad;
        std::vector< std::vector<float> > nlgrad_upd;
        
        {
            std::vector<float> aux;
            for(int i = 0; i < dim; i++)
            {
                aux.assign((int) wxy[i].size(), 0.0f);
                
                q.push_back(aux);
                q_upd.push_back(aux);
                nlgrad.push_back(aux);
                nlgrad_upd.push_back(aux);
            }
        }
        
        float *nldiv = new float[dim];
        float *nldiv_upd = new float[dim];
        
        fpClear(nldiv, 0.0f, dim);
        fpClear(nldiv_upd, 0.0f, dim);
        
        // Adaptive primal-dual Chambolle-Pock algorithm
        float alpha = 0.5f;
        float beta = 0.95f;
        float delta = 1.5f;
        float eta = 0.95f;
        float gamma = 0.75f;
        float s_rang = 255.0f;
        tau = 0.5f;
        sigma = 0.5f;
        
        float alpha0 = alpha;
        float tolPrimal = tol * dim;
        float tolDual = tolPrimal * 2.0f; // take into account number of dual variables
        float primalRes = tolPrimal;
        float dualRes = tolDual;
        int iter = 0;
        
        while(iter < maxIter && primalRes >= tolPrimal && dualRes >= tolDual)
        {
            // Compute proximity operator of primal energy
            proxGL2(u_upd, u, f, mask, 1.0f, tau, div, nldiv, dim);
            
            // Compute gradient, proximity operator and divergence for TV regularization term
            forward_gradient(u_upd, gradx_upd, grady_upd, width, height);
            proxTV(px_upd, py_upd, px, py, gradx_upd, grady_upd, gradx, grady, lmbTV, sigma, dim);
            divergence(px_upd, py_upd, div_upd, width, height);
            
            // Compute gradient, proximity operator and divergence for NL regularization term
            nl_gradient(u_upd, nlgrad_upd, wxy, posxy, dim);
            proxNL(q_upd, q, nlgrad_upd, nlgrad, lmbNL, sigma, dim);
            nl_divergence(q_upd, nldiv_upd, wxy, wyx, posyx, posw, dim);
            
            // Backtracking adaptive PDHG
            primalRes = 0.0f;
            dualRes = 0.0f;
            float unorm = 0.0f;
            float pnorm = 0.0f;
            float qnorm = 0.0f;
            float inprod = 0.0f;
            
            // Compute primal and dual residuals
            for(int i = 0; i < dim; i++)
            {
                float udif = u[i] - u_upd[i];
                float pxdif = px[i] - px_upd[i];
                float pydif = py[i] - py_upd[i];
                
                float gradxdif = gradx[i] - gradx_upd[i];
                float gradydif = grady[i] - grady_upd[i];
                
                float divdif = div[i] - div_upd[i];
                float nldivdif = nldiv[i] - nldiv_upd[i];
                
                unorm += (udif * udif);
                pnorm += (pxdif * pxdif + pydif * pydif);
                inprod += (gradxdif * pxdif + gradydif * pydif);
                
                primalRes += fabs((udif / tau) + divdif + nldivdif);
                
                float xval = (pxdif / sigma) - gradxdif;
                float yval = (pydif / sigma) - gradydif;
                dualRes += sqrt(xval * xval + yval * yval);
                
                float qvalsum = 0.0f;
                
                for(int w = 0; w < (int) wxy[i].size(); w++)
                {
                    float qdif = q[i][w] - q_upd[i][w];
                    float nlgraddif = nlgrad[i][w] - nlgrad_upd[i][w];
                    
                    qnorm += (qdif * qdif);
                    inprod += (nlgraddif * qdif);
                    
                    float qval = (qdif / sigma) - nlgraddif;
                    qvalsum = qval * qval;
                }
                dualRes += sqrt(qvalsum);

            }
            
            
            // Backtracking condition
            float b = (2.0f * tau * sigma * inprod) / (gamma * sigma * unorm + gamma * tau * (pnorm + qnorm));
            
            if(b > 1) // Backtracking condition doesn't hold, thus reduce both step sizes
            {
                tau = (beta * tau) / b;
                sigma = (beta * sigma) / b;
                alpha = alpha0;
                
                //if(flagVerbose)
                    printf("... Backtracking condition activated at iteration %i ...\n", iter);
                
            } else
            {
                float dualResInf = s_rang * delta * dualRes;
                float dualResSup = (s_rang * dualRes) / delta;
                
                if(primalRes > dualResInf) // Primal residual larger than dual residual, thus increase primal stepsize
                {
                    tau = tau / (1.0f - alpha);
                    sigma = sigma * (1.0f - alpha);
                    alpha = alpha * eta;
                    
                } else if(primalRes < dualResSup) // Dual residual larger than primal residual, thus increase dual stepsize
                {
                    tau = tau * (1.0f - alpha);
                    sigma = sigma / (1.0f - alpha);
                    alpha = alpha * eta;
                }
                
                // Update variables
                fpCopy(u_upd, u, dim);
                fpCopy(px_upd, px, dim);
                fpCopy(py_upd, py, dim);
                fpCopy(gradx_upd, gradx, dim);
                fpCopy(grady_upd, grady, dim);
                fpCopy(div_upd, div, dim);
                fpCopy(nldiv_upd, nldiv, dim);
                
                for(int i = 0; i < dim; i++)
                {
                    q[i].assign(q_upd[i].begin(), q_upd[i].end());
                    nlgrad[i].assign(nlgrad_upd[i].begin(), nlgrad_upd[i].end());
                }
                
                iter++;
                
                //if(flagVerbose)
                    printf("Iter %5.0i - PrimalRes = %f, DualRes = %f, tau = %f, sigma = %f\n", iter, primalRes, dualRes, tau, sigma);
            }
        }
        
        // Delete allocated memory
        delete[] u_upd;
        delete[] px; delete[] py; delete[] px_upd; delete[] py_upd;
        delete[] gradx; delete[] grady; delete[] gradx_upd; delete[] grady_upd;
        delete[] div; delete[] div_upd; delete[] nldiv; delete[] nldiv_upd;
    }
    

    
    // Compute proximity operator of L2 data-fitting term
    void proxGL2(float *u_upd, float *u, float *f, float *mask, float lmbFit, float tau,
               float *div, float *nldiv, int dim)
    {
        for(int i = 0; i < dim; i++)
            u_upd[i] = (u[i] + tau * (div[i] + nldiv[i] + lmbFit * mask[i] * f[i])) / (1.0f + lmbFit * tau * mask[i]);
    }
    
    
    
    // Compute proximity operator of L1 data-fitting term
    void proxGL1(float *u_upd, float *u, float *f, float *mask, float lmbFit, float tau,
                 float *div, float *nldiv, int dim)
    {
        for(int i = 0; i < dim; i++)
        {
            float arg = u[i] + tau * (div[i] + nldiv[i]);
            float thres = tau * lmbFit * mask[i];
            
            if(arg < f[i] - thres)
            {
                u_upd[i] = arg + thres;
                
            } else if(arg > f[i] + thres)
            {
                u_upd[i] = arg - thres;
                
            } else
            {
                u_upd[i] = f[i];
            }
        }
    }
    
    
    
    // Compute proximity operator of TV regularization
    void proxTV(float *px_upd, float *py_upd, float *px, float *py, float *gradx_upd, float *grady_upd,
                float *gradx, float *grady, float lmbTV, float sigma, int dim)
    {
        
        for(int i = 0; i < dim; i++)
        {
            float argx = px[i] + sigma * (2.0f * gradx_upd[i] - gradx[i]);
            float argy = py[i] + sigma * (2.0f * grady_upd[i] - grady[i]);
            
            float norm = sqrt(argx * argx + argy * argy);
            float maxval = MAX(lmbTV, norm);
            
            px_upd[i] = (lmbTV * argx) / maxval;
            py_upd[i] = (lmbTV * argy) / maxval;
        }
    }
    
    
    
    // Compute proximity operator of NL regularization
    void proxNL(std::vector< std::vector <float> > &q_upd, std::vector< std::vector <float> > &q,
                std::vector< std::vector <float> > &nlgrad_upd, std::vector< std::vector <float> > &nlgrad,
                float lmbNL, float sigma, int dim)
    {
        
        for(int i = 0; i < dim; i++)
        {
            float norm = 0.0f;
            
            for(int w = 0; w < (int) q[i].size(); w++)
            {
                float nlarg = q[i][w] + sigma * (2.0f * nlgrad_upd[i][w] - nlgrad[i][w]);
                norm += (nlarg * nlarg);
                
                q_upd[i][w] = lmbNL * nlarg;
            }
            
            float maxval = MAX(lmbNL, sqrt(norm));
            
            for(int w = 0; w < (int) q[i].size(); w++)
                q_upd[i][w] /= maxval;
        }
    }
    
    
    
    // Compute local gradient at each pixel via forward differences
    void forward_gradient(float *f, float *fx, float *fy, int nx, int ny)
    {
        for(int i = 0; i < ny-1; i++)
        {
            int p, p1, p2;
            
            for(int j = 0; j < nx-1; j++)
            {
                p  = i * nx + j;
                p1 = p + 1;
                p2 = p + nx;
                fx[p] = f[p1] - f[p];
                fy[p] = f[p2] - f[p];
            }
        }
        
        
        // Compute the gradient on the last row
        for(int j = 0; j < nx-1; j++)
        {
            const int p = (ny-1) * nx + j;
            fx[p] = f[p+1] - f[p];
            fy[p] = 0;
        }
        
        
        // Compute the gradient on the last column
        for(int i = 1; i < ny; i++)
        {
            const int p = i * nx-1;
            fx[p] = 0;
            fy[p] = f[p+nx] - f[p];
        }
        
        fx[ny * nx - 1] = 0;
        fy[ny * nx - 1] = 0;
    }
    
    
    
    // Compute local divergence at each pixel as minus the adjoint operator of the local gradient, div = -nabla^T
    void divergence(float *v1, float *v2, float *div, int nx, int ny)
    {
        for(int i = 1; i < ny-1; i++)
            for(int j = 1; j < nx-1; j++)
            {
                const int p  = i * nx + j;
                const int p1 = p - 1;
                const int p2 = p - nx;
                
                const double v1x = v1[p] - v1[p1];
                const double v2y = v2[p] - v2[p2];
                
                div[p] = v1x + v2y;
            }
        
        // Compute the divergence on the first and last rows
        for (int j = 1; j < nx-1; j++)
        {
            const int p = (ny-1) * nx + j;
            div[j] = v1[j] - v1[j-1] + v2[j];
            div[p] = v1[p] - v1[p-1] - v2[p-nx];
        }
        
        
        // Compute the divergence on the first and last columns
        for (int i = 1; i < ny-1; i++)
        {
            const int p1 = i * nx;
            const int p2 = (i+1) * nx - 1;
            div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
            div[p2] = -v1[p2-1] + v2[p2] - v2[p2 - nx];
        }
        
        div[0]         =  v1[0]         + v2[0];
        div[nx-1]      = -v1[nx - 2]    + v2[nx - 1];
        div[(ny-1)*nx] =  v1[(ny-1)*nx] - v2[(ny-2)*nx];
        div[ny*nx-1]   = -v1[ny*nx - 2] - v2[(ny-1)*nx - 1];
    }
    
    
    
    // Compute nonlocal gradient operator at each pixel
    void nl_gradient(float *data, std::vector< std::vector<float> > &nlgrad, std::vector< std::vector<float> > &weights,
                     std::vector< std::vector<int> > &positions, int dim)
    {
        for(int l = 0; l < dim; l++)
        {
            for(int w = 0; w < (int) weights[l].size(); w++)
            {
                int l0 = positions[l][w];
                nlgrad[l][w] = (data[l0] - data[l]) * weights[l][w];
            }
        }
    }
    
    
    
    // Compute nonlocal divergence at each pixel as minus the adjoint operator of the nonlocal gradient, div_w = -nabla_w^T
    void nl_divergence(std::vector< std::vector<float> > &data, float *nldiv, std::vector< std::vector<float> > &wxy,
                       std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx,
                       std::vector< std::vector<int> > &posw, int dim)
    {
        for(int l = 0; l < dim; l++)
        {
            float divl = 0.0f;
            
            // Terms involving w(x,y)
            for(int w = 0; w < (int) wxy[l].size(); w++)
                divl += data[l][w] * wxy[l][w];
            
            // Terms involving w(y,x)
            for(int w = 0; w < (int) wyx[l].size(); w++)
            {
                int l0 = posyx[l][w];
                int w0 = posw[l][w];
                divl -= data[l0][w0] * wyx[l][w];
            }
            
            nldiv[l] = divl;
        }
    }
    
    
    
    // Compute error between two consecutive iterations
    float relative_error(float *u, float *u_upd, int dim)
    {
        float error = 0.0f;
        
        for(int i = 0; i < dim; i++)
        {
            float value = u_upd[i] - u[i];
            error += (value * value);
        }
        
        error /= ((float) dim);
        error = sqrtf(error);
        
        return error;
    }
    
    
    
    
    // Compute weights for NL regularization on a prescribed image using all neighbouring pixels
    void nlweights_reg(float **image, std::vector< std::vector<float> > &wxy, std::vector< std::vector<int> > &posxy,
                       std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx,
                       std::vector< std::vector<int> > &posw, float hSim, float hClose, int numNeighbours,
                       int reswind, int compwind, int flagNormalize, int num_channels, int width, int height)
    {
        
        // Research window size
        int resdim = (2 * reswind + 1) * (2 * reswind + 1);
        
        // Tabulate the function Exp(-x) for x>0.
        int luttaille = (int) (LUTMAX * LUTPRECISION);
        float *lut = new float[luttaille];
        libUSTG::wxFillExpLut(lut, luttaille);
        
        // Auxiliar vectors
        float *dist = new float[resdim];
        float *pos = new float[resdim];
        
        // For each pixel l, compute weights w(l,l0)
        //#pragma omp parallel shared(image, lut, wxy, posxy, numxy)
        {
            //#pragma omp for schedule(dynamic) nowait
            for(int y = 0; y < height; y++)
            {
                // Fix row
                int ly = y * width;
                
                for(int x = 0; x < width; x++)
                {
                    // Index of current pixel
                    int l = ly + x;
                    
                    // Reduce window size if we are close to boundary
                    int compwind0 = MAX(MIN(compwind, MIN(width-1-x-reswind, MIN(height-1-y-reswind, MIN(x-reswind, y-reswind)))), 0);
                    
                    // Learning zone depending on the window size
                    int imin = MAX(x - reswind, compwind0);
                    int jmin = MAX(y - reswind, compwind0);
                    int imax = MIN(x + reswind, width - compwind0 - 1);
                    int jmax = MIN(y + reswind, height - compwind0 - 1);
                    
                    // Adapt filter parameters to size of comparison window
                    float compdim = (float) (2 * compwind0 + 1) * (2 * compwind0 + 1);
                    float hSim2 = hSim * hSim * compdim;
                    float hClose2 = hClose * hClose;
                    
                    // Clean vectors
                    libUSTG::fpClear(dist, 0.0f, resdim);
                    libUSTG::fpClear(pos, 0.0f, resdim);
                    
                    // Auxiliar variables
                    float dMin = fLarge;
                    int wcentral = 0;
                    int w = 0;
                    
                    // Compute distance for each pixel in the resdim-neighbourhood
                    for(int j = jmin; j <= jmax; j++)
                    {
                        int l0y = j * width;
                        
                        for(int i = imin; i <= imax; i++)
                        {
                            // Compute patch-based similarity distance
                            float dSim = libUSTG::fiL2FloatDist(image, image, x, y, i, j, compwind0, compwind0, num_channels,
                                                                width, width);
                            dSim /= hSim2;
                            
                            // Compute spatial distance
                            float dClose = (float) ((x - i) * (x - i) + (y - j) * (y - j));
                            dClose /= hClose2;
                            
                            // Save bilateral distance
                            dist[w] = dSim + dClose;
                            
                            // Save w-position of central pixel in neighbourhood and the minimum bilateral distance
                            if((i == x) && (j == y))
                            {
                                wcentral = w;
                                
                            } else
                            {
                                if(dist[w] < dMin)
                                    dMin = dist[w];
                            }
                            
                            // Save neighbouring pixel position
                            pos[w] = (float) (l0y + i);
                            
                            // Update index
                            w++;
                        }
                    }
                    
                    // Assign minimum distance to central pixel
                    dist[wcentral] = dMin;
                    
                    // Adapt number of neighbouring pixels to window size
                    int numNeighbours0 = MIN(numNeighbours, w);
                    
                    // Order bilateral distances
                    libUSTG::fpQuickSort(dist, pos, w);
                    
                    // Auxiliar vectors
                    std::vector<float> auxw;
                    std::vector<int> auxp;
                    
                    // Compute weight w(l,l0) for each l0
                    float weight_sum = 0.0f;
                    
                    for(int r = 0; r < numNeighbours0; r++)
                    {
                        float weight = libUSTG::wxSLUT(dist[r], lut);
                        auxw.push_back(weight);
                        weight_sum += weight;
                        
                        auxp.push_back((int) pos[r]);
                    }
                    
                    // Normalize weights w(l,l0) with respect to l0
                    if(flagNormalize)
                    {
                        if(weight_sum > fTiny)
                        {
                            for(int r = 0; r < numNeighbours0; r++)
                            {
                                float wval = auxw[r];
                                auxw[r] = sqrt(wval / weight_sum);
                            }
                        } else
                        {
                            auxw.clear();
                            auxp.clear();
                            
                            auxw.push_back(1.0f);
                            auxp.push_back(l);
                        }
                        
                    } else
                    {
                        for(int r = 0; r < numNeighbours0; r++)
                        {
                            float wval = auxw[r];
                            auxw[r] = sqrt(wval);
                        }
                    }
                    
                    wxy.push_back(auxw);
                    posxy.push_back(auxp);
                }
            }
        }
        
        // For each pixel l, compute weights w(l0,l)
        {
            //#pragma omp for schedule(dynamic) nowait
            for(int y = 0; y < height; y++)
            {
                int ly = y * width;
                
                for(int x = 0; x < width; x++)
                {
                    // Current central pixel
                    int l = ly + x;
                    
                    // Reduce window size if we are close to boundary
                    int compwind0 = MAX(MIN(compwind, MIN(width-1-x-reswind, MIN(height-1-y-reswind, MIN(x-reswind, y-reswind)))), 0);
                    
                    // Learning zone depending on the window size
                    int imin = MAX(x - reswind, compwind0);
                    int jmin = MAX(y - reswind, compwind0);
                    int imax = MIN(x + reswind, width - compwind0 - 1);
                    int jmax = MIN(y + reswind, height - compwind0 - 1);
                    
                    // Auxiliar vectors
                    std::vector<float> auxw;
                    std::vector<int> auxp;
                    std::vector<int> auxpw;
                    
                    // Look at research window around l
                    for(int j = jmin; j <= jmax; j++)
                    {
                        int l0j = j * width;
                        
                        for(int i = imin; i <= imax; i++)
                        {
                            // Neighbouring pixel
                            int l0 = l0j + i;
                            
                            // Consider pixels l0 such that w(l0,l)>0
                            for(int w = 0; w < (int) wxy[l0].size(); w++)
                            {
                                if(posxy[l0][w] == l)
                                {
                                    float weight = wxy[l0][w];
                                    
                                    if(weight > fTiny)
                                    {
                                        auxw.push_back(weight);
                                        auxp.push_back(l0);
                                        auxpw.push_back(w);
                                    }
                                    
                                    break;
                                }
                            }
                        }
                    }
                    
                    // Save outputs
                    wyx.push_back(auxw);
                    posyx.push_back(auxp);
                    posw.push_back(auxpw);
                }
            }
        }
        
        // Delete allocated memory
        delete[] lut;
        delete[] dist;
        delete[] pos;
    }
    
    
    
    // Compute weights for NL regularization on a prescribed image using only neighbouring pixels activated by a mask
    void nlweights_reg_mask(float **image, float *mask, std::vector< std::vector<float> > &wxy,
                            std::vector< std::vector<int> > &posxy, std::vector< std::vector<float> > &wyx,
                            std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw,
                            float hSim, float hClose, int numNeighbours, int reswind, int compwind, int flagNormalize,
                            int num_channels, int width, int height)
    {
        
        // Research window size
        int resdim = (2 * reswind + 1) * (2 * reswind + 1);
        
        // Tabulate the function Exp(-x) for x>0.
        int luttaille = (int) (LUTMAX * LUTPRECISION);
        float *lut = new float[luttaille];
        libUSTG::wxFillExpLut(lut, luttaille);
        
        // Auxiliar vectors
        float *dist = new float[resdim];
        float *pos = new float[resdim];
        
        // For each pixel l, compute weights w(l,l0)
        //#pragma omp parallel shared(image, lut, wxy, posxy, numxy)
        {
            //#pragma omp for schedule(dynamic) nowait
            for(int y = 0; y < height; y++)
            {
                // Fix row
                int ly = y * width;
                
                for(int x = 0; x < width; x++)
                {
                    // Index of current pixel
                    int l = ly + x;
                    
                    // Reduce window size if we are close to boundary
                    int compwind0 = MAX(MIN(compwind, MIN(width-1-x-reswind, MIN(height-1-y-reswind, MIN(x-reswind, y-reswind)))), 0);
                    
                    // Learning zone depending on the window size
                    int imin = MAX(x - reswind, compwind0);
                    int jmin = MAX(y - reswind, compwind0);
                    int imax = MIN(x + reswind, width - compwind0 - 1);
                    int jmax = MIN(y + reswind, height - compwind0 - 1);
                    
                    // Adapt filter parameters to size of comparison window
                    float compdim = (float) (2 * compwind0 + 1) * (2 * compwind0 + 1);
                    float hSim2 = hSim * hSim * compdim;
                    float hClose2 = hClose * hClose;
                    
                    // Auxiliar variables
                    float dMin = fLarge;
                    int wcentral = 0;
                    int w = 0;
                    
                    // Clean vectors
                    libUSTG::fpClear(dist, 0.0f, resdim);
                    libUSTG::fpClear(pos, 0.0f, resdim);
                    
                    // Compute distance for each pixel in the resdim-neighbourhood
                    for(int j = jmin; j <= jmax; j++)
                    {
                        int l0y = j * width;
                        
                        for(int i = imin; i <= imax; i++)
                        {
                            int l0 = j * width + i;
                            
                            // Compare only with neighbouring pixels activated by the mask
                            if(mask[l0] == 1)
                            {
                                // Compute patch-based similarity distance
                                float dSim = libUSTG::fiL2FloatDist(image, image, x, y, i, j, compwind0, compwind0, num_channels,
                                                                    width, width);
                                dSim /= hSim2;
                                
                                // Compute spatial distance
                                float dClose = (float) ((x - i) * (x - i) + (y - j) * (y - j));
                                dClose /= hClose2;
                                
                                // Save bilateral distance
                                dist[w] = dSim + dClose;
                                
                                // Save w-position of central pixel in neighbourhood and the minimum bilateral distance
                                if((i == x) && (j == y))
                                {
                                    wcentral = w;
                                    
                                } else
                                {
                                    if(dist[w] < dMin)
                                        dMin = dist[w];
                                }
                                
                                // Save neighbouring pixel position
                                pos[w] = (float) (l0y + i);
                                
                                // Update index
                                w++;
                            }
                        }
                    }
                    
                    // Assign minimum distance to central pixel
                    if(mask[l] == 1)
                        dist[wcentral] = dMin;
                    
                    // Adapt number of neighbouring pixels to window size
                    int numNeighbours0 = MIN(numNeighbours, w);
                    
                    // Order bilateral distances
                    libUSTG::fpQuickSort(dist, pos, w);
                    
                    // Auxiliar vectors
                    std::vector<float> auxw;
                    std::vector<int> auxp;
                    
                    // Compute weight w(l,l0) for each l0
                    float weight_sum = 0.0f;
                    
                    for(int r = 0; r < numNeighbours0; r++)
                    {
                        float weight = libUSTG::wxSLUT(dist[r], lut);
                        auxw.push_back(weight);
                        weight_sum += weight;
                        
                        auxp.push_back((int) pos[r]);
                    }
                    
                    // Normalize weights w(l,l0) with respect to l0
                    if(flagNormalize)
                    {
                        if(weight_sum > fTiny)
                        {
                            for(int r = 0; r < numNeighbours0; r++)
                            {
                                float wval = auxw[r];
                                auxw[r] = sqrt(wval / weight_sum);
                            }
                        } else
                        {
                            auxw.clear();
                            auxp.clear();
                            
                            auxw.push_back(1.0f);
                            auxp.push_back(l);
                        }
                        
                    } else
                    {
                        for(int r = 0; r < numNeighbours0; r++)
                        {
                            float wval = auxw[r];
                            auxw[r] = sqrt(wval);
                        }
                    }
                    
                    wxy.push_back(auxw);
                    posxy.push_back(auxp);
                }
            }
        }
        
        // For each pixel l, compute weights w(l0,l)
        {
            //#pragma omp for schedule(dynamic) nowait
            for(int y = 0; y < height; y++)
            {
                int ly = y * width;
                
                for(int x = 0; x < width; x++)
                {
                    // Current central pixel
                    int l = ly + x;
                    
                    // Reduce window size if we are close to boundary
                    int compwind0 = MAX(MIN(compwind, MIN(width-1-x-reswind, MIN(height-1-y-reswind, MIN(x-reswind, y-reswind)))), 0);
                    
                    // Learning zone depending on the window size
                    int imin = MAX(x - reswind, compwind0);
                    int jmin = MAX(y - reswind, compwind0);
                    int imax = MIN(x + reswind, width - compwind0 - 1);
                    int jmax = MIN(y + reswind, height - compwind0 - 1);
                    
                    // Auxiliar vectors
                    std::vector<float> auxw;
                    std::vector<int> auxp;
                    std::vector<int> auxpw;
                    
                    // Look at research window around l
                    for(int j = jmin; j <= jmax; j++)
                    {
                        int l0j = j * width;
                        
                        for(int i = imin; i <= imax; i++)
                        {
                            // Neighbouring pixel
                            int l0 = l0j + i;
                            
                            // Consider pixels l0 such that w(l0,l)>0
                            for(int w = 0; w < (int) wxy[l0].size(); w++)
                            {
                                if(posxy[l0][w] == l)
                                {
                                    float weight = wxy[l0][w];
                                    
                                    if(weight > fTiny)
                                    {
                                        auxw.push_back(weight);
                                        auxp.push_back(l0);
                                        auxpw.push_back(w);
                                    }
                                    
                                    break;
                                }
                            }
                        }
                    }
                    
                    // Save outputs
                    wyx.push_back(auxw);
                    posyx.push_back(auxp);
                    posw.push_back(auxpw);
                }
            }
        }
        
        // Delete allocated memory
        delete[] lut;
        delete[] dist;
        delete[] pos;
    }

} // libUSTGFiltVar


