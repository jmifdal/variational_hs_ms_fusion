#ifndef _LIBHYPERSPECTRALWEIGHTS_H_
#define _LIBHYPERSPECTRALWEIGHTS_H_

#include "../library/libBasic.h"
#include "../library/libImage.h"

//#include "libBasic.h"
//#include "libImage.h"

using namespace std;

namespace libUSTGHYPER
{
    
    // Values and positions
    struct ValAndPos {
        float value;
        int position;
    };
    
    
    // Values, positions and comparison values
    struct ValPosCompVect {
        float value;
        float comp_value;
        int position;
    };
    
    // Computes the sqrt(weight) on the MS data, so 3D patches are considered
    // The patch distance across channels is coupled with the L2 norm.
    // The same weight distribution is obtained for all HS channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    void nlweights_ms(float **mData, std::vector< std::vector<float> > &wxy, std::vector< std::vector<int> > &posxy,
                      std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw,
                      float hPatch, float hSpatial, int numNeighbours, int reswind, int compwind, int flagNormalize, int ms_channels,
                      int width, int height);
    
    // Computes the sqrt(weight) on the MS data, so 3D patches are considered.
    // The contribution of each MS channel to the patch-based distance is weighted by the coefficients of the spectral downsampling matrix.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    void nlweights_ms_Smatrix_new(libUSTG::cflimage ms_image,float **Smatrix ,
    		std::vector< std::vector<float> > &wxy,std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
    		std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, int hs_channels, float hSim,
    		float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize);

    // Computes the sqrt(weight) on the MS data, so 3D patches are considered.
    // The contribution of each MS channel to the patch-based distance is weighted by the coefficients of the spectral downsampling matrix.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    void nlweights_ms_Smatrix(float **mData, float **Smatrix, std::vector< std::vector< std::vector<float> > > &wxy,
                              std::vector< std::vector< std::vector<int> > > &posxy, std::vector< std::vector< std::vector<float> > > &wyx,
                              std::vector< std::vector< std::vector<int> > > &posyx, std::vector< std::vector< std::vector<int> > > &posw,
                              float hPatch, float hSpatial, int numNeighbours, int reswind, int compwind, int flagNormalize,
                              int hs_channels, int ms_channels, int width, int height);
    
    // For each HS channel, computes the sqrt(weight) on the corresponding upsampled HS image, so 2D patches are considered.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
    void nlweights_hs_band(float **hDataInt, std::vector< std::vector< std::vector<float> > > &wxy,
                           std::vector< std::vector< std::vector<int> > > &posxy, std::vector< std::vector< std::vector<float> > > &wyx,
                           std::vector< std::vector< std::vector<int> > > &posyx, std::vector< std::vector< std::vector<int> > > &posw,
                           float hPatch, float hSpatial, int numNeighbours, int reswind, int compwind, int flagNormalize, int hs_channels,
                           int width, int height);
    
    // Computes the sqrt(weights) by first choosing most similar neighbours on the full MS image (3D patches) and computing then
    // the weights using 2D patches centered at the pre-selected pixels on the corresponding band.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h*dim+i][n] for each channel h, each pixel i, and each neighbour n.
    void nlweights_ms_3Dclassify_2Dband(libUSTG::cflimage hyperDataUpsampled, libUSTG::cflimage ms_image,std::vector< std::vector<float> > &wxy,
    		std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
    		std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, float hSim,
    		float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize);


    //In this function we classify 2D patch of the hyperspectral image based 3D
    //patches computed on the multispectral image and weighted with the Smatrix
    
    //TODO : Joan, maybe you should add tests to verify that the upsampled hyperspectral image and the multispectral image
    //have the same spatial resolution. And also test if the number of the hyperspectral bands of the upsampled hyperspectral image and
    //those of the multispectral bands are consistent of the dimensions of the Smatrix
    void nlweights_ms_3Dclassify_Smatrix_2Dband(libUSTG::cflimage hyperDataUpsampled,libUSTG::cflimage ms_image,float **Smatrix ,
                                                std::vector< std::vector<float> > &wxy,std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
                                                std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, float hSim,
                                                float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize);
    
    //This routines computes 3D distances on the multispectral image weighted with coefficient of the Smatrix
    //the distances are unclassified
    void nlweights_3Dall_ms_Smatrix_unclassified(float *ms_image_new_ptr, std::vector<float> &vect_xy, float hSim, float hClose,
                                                 int reswind, int spat_patch_step, int ms_channels,int width_ms_new,
                                                 int height_ms_new,float **Smatrix,int x, int y, int num_hs_channel);
    
    // Computes the sqrt(weights) on the full upsampled HS image, so 3D patches are considered.
    // The patch distance across channels is coupled with the L2 norm.
    // The same weight distribution is obtained for all HS channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    void nlweights_hs_3Dall(libUSTG::cflimage hyperDataUpsampled, std::vector< std::vector<float> > &wxy,
                            std::vector< std::vector<int> > &posxy, std::vector< std::vector<float> > &wyx,
                            std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw, float hSim, float hClose,
                            int numSimPixels, int reswind, int spat_patch_step, int flagNormalize);
    
    // Computes the sqrt(weights) using the nearest bands w.r.t. the central band, so 3D patches are considered.
    // The patch distance across channels is coupled with the L2 norm.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h*dim+i][n] for each channel h, each pixel i, and each neighbour n.
    void nl_weights_hs_3Dnearest(libUSTG::cflimage hyperDataUpsampled, std::vector< std::vector<float> > &wxy,
                                 std::vector< std::vector<int> > &posxy, std::vector< std::vector<float> > &wyx,
                                 std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw, float hSim, float hClose,
                                 int numSimPixels, int reswind, int spat_patch_step, int spec_patch_step, int flagNormalize);
    
    // Computes the sqrt(weights) by first choosing most similar neighbours on the full upsampled HS image (3D patches) and computing then
    // the weights using 2D patches centered at the pre-selected pixels on the corresponding band.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h*dim+i][n] for each channel h, each pixel i, and each neighbour n.
    void nlweights_hs_3Dclassify_2Dband(libUSTG::cflimage hyperDataUpsampled, std::vector< std::vector<float> > &wxy,
                                        std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
                                        std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, float hSim,
                                        float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize);
    
    // Computes weights on the full upsampled HS image, so 3D patches are considered, but without classifiying them.
    // The same weight is obtained for all HS channels.
    // Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
    void nlweights_3Dall_unclassified(float *im_3d_new_ptr, std::vector< std::vector<float> > &wxy, float hSim, float hClose,int reswind,
                                      int spat_patch_step, int hs_channels, int width_new, int height_new);
    
    // Sorts two vectors wrt the values vector (ascending order)
    void sort_vects(std::vector<float> &values_vect, std::vector<int> &positions_vect);
    
    // Sort with respect to a comparison vector (ascending order)
    void sort_wrt_comp_vect(std::vector<float> &values_vect, std::vector<int> &positions_vect, const std::vector<float> &comparison_vect);
    
    void array_slicing_coef(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height, int width,
                            float *im_to_slice, int ind, std::vector<float>  &sliced_vector,float coef);
    
    void array_slicing(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height, int width,
                       float *im_to_slice, std::vector<float>  &sliced_vector);

    float sum_values_vector(std::vector<float> &vec, int size_vec);
    
    void substract_two_vectors(std::vector<float> &vec1, std::vector<float> &vec2, int size_vec, std::vector<float> &res_vector);
    
    void power_two_vector(std::vector<float> &vector, int size_vec, std::vector<float> &res_vector);
    
    // Prints the values of a 1D pointer array.
    void print_values_array_1D_pt(float *array, int size);
    
    // Prints the values of an array from the cflimage class.
    void print_values_array(libUSTG::cflimage vec);
    
    // Extends the SPATIAL dimensions of an image by mirror effect which helps in case the treatment of the image takes the borders into
    // account. The image to extend should be a 3D image.
    libUSTG::cflimage extend_image_dimensions_spat(libUSTG::cflimage im_to_extend, int spat_patch_step);
    
    // Extends the SPECTRAL dimensions of an image by mirror effect which helps in case the treatment of the image takes the borders into
    // account. The image to extend should be a 3d image.
    libUSTG::cflimage  extend_image_dimensions_spec(libUSTG::cflimage im_to_extend, int spec_patch_step);
    
} // libUSTGHYPER


#endif
