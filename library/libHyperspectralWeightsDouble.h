#ifndef _LIBHYPERSPECTRALWEIGHTSDOUBLE_H_
#define _LIBHYPERSPECTRALWEIGHTSDOUBLE_H_

#include "../library/libBasic.h"
#include "../library/libImage.h"


using namespace std;

namespace libUSTGHYPERDOUBLE
{
    
    // Values and positions
    struct ValAndPos {
        double value;
        int position;
    };
    
    
    // Values, positions and comparison values
    struct ValPosCompVect {
        double value;
        double comp_value;
        int position;
    };
    
    // Classify 2D patch of the hyperspectral image based 3D patches computed on the multispectral image and weighted with the Smatrix
    // TODO : Joan, maybe you should add tests to verify that the upsampled hyperspectral image and the multispectral image
    // have the same spatial resolution. And also test if the number of the hyperspectral bands of the upsampled hyperspectral image and
    // those of the multispectral bands are consistent of the dimensions of the Smatrix
    void nlweights_ms_3Dclassify_Smatrix_2Dband_double(libUSTG::cflimage hyperDataUpsampled, libUSTG::cflimage ms_image, double **Smatrix ,
                                                       std::vector< std::vector<double> > &wxy,std::vector< std::vector<int> > &pxy,
                                                       std::vector< std::vector<double> > &wyx, std::vector< std::vector<int> > &pyx,
                                                       std::vector< std::vector<int> > &posw, double hSim,  double hClose, int numSimPixels,
                                                       int reswind, int spat_patch_step, int flagNormalize);
    
    // Compute 3D distances on the multispectral image weighted with coefficient of the Smatrix.
    // The distances are unclassified.
    // Width_ms_new and height_ms_new are the spatial dimensions of the multispectral image after spatial extension.
    void nlweights_3Dall_ms_Smatrix_unclassified_double(double *ms_image_new_ptr, std::vector<double> &vect_xy, double hSim, double hClose,
                                                        int reswind, int spat_patch_step, int ms_channels, int width_ms_new, int height_ms_new,
                                                        double **Smatrix, int x, int y, int num_hs_channel);
    
    // Sort with respect to a comparison vector (ascending order)
    void sort_wrt_comp_vect_double(std::vector<double> &values_vect, std::vector<int> &positions_vect, const std::vector<double> &comparison_vect);
    
    void array_slicing_coef_double(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height,
                                   int width, double *im_to_slice, int ind, std::vector<double> &sliced_vector, double coef);
    
    void array_slicing_double(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height,
                              int width, double *im_to_slice, std::vector<double> &sliced_vector);
    
    double sum_values_vector_double(std::vector<double> &vec, int size_vec);
    
    void substract_two_vectors_double(std::vector<double> &vec1, std::vector<double> &vec2, int size_vec, std::vector<double> &res_vector);
    
    void power_two_vector_double(std::vector<double> &vector, int size_vec, std::vector<double> &res_vector);
    
    libUSTG::cflimage extend_image_dimensions_spat(libUSTG::cflimage im_to_extend, int spat_patch_step);
    
} // libUSTGHYPERDOUBLE


#endif
