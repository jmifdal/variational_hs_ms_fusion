#include "libHyperspectralWeightsDouble.h"


namespace libUSTGHYPERDOUBLE
{

    // Classify 2D patch of the hyperspectral image based 3D patches computed on the multispectral image and weighted with the Smatrix
    // TODO : Joan, maybe you should add tests to verify that the upsampled hyperspectral image and the multispectral image
    // have the same spatial resolution. And also test if the number of the hyperspectral bands of the upsampled hyperspectral image and
    // those of the multispectral bands are consistent of the dimensions of the Smatrix
    void nlweights_ms_3Dclassify_Smatrix_2Dband_double(libUSTG::cflimage hyperDataUpsampled, libUSTG::cflimage ms_image, double **Smatrix ,
                                                       std::vector< std::vector<double> > &wxy,std::vector< std::vector<int> > &pxy,
                                                       std::vector< std::vector<double> > &wyx, std::vector< std::vector<int> > &pyx,
                                                       std::vector< std::vector<int> > &posw, double hSim,  double hClose, int numSimPixels,
                                                       int reswind, int spat_patch_step, int flagNormalize)
    {
        
        // Tabulate the function Exp(-x) for x>0.
        int luttaille = (int) (dLUTMAX * dLUTPRECISION);
        double *lut = new double[luttaille];
        libUSTG::wxFillExpLut(lut, luttaille);
        
        // Extension of the images. In this case the image should be extended only spatially
        libUSTG::cflimage im_3d_new = extend_image_dimensions_spat(hyperDataUpsampled, spat_patch_step); // spatial extension
        int size = im_3d_new.whc();
        double *im_3d_new_ptr = new double[size];
 
        for(int i = 0; i < size; i++)
            im_3d_new_ptr[i] = (double) (im_3d_new.v()[i]);
        
        libUSTG::cflimage ms_3d_new = extend_image_dimensions_spat(ms_image, spat_patch_step); // spatial extension
        size = ms_3d_new.whc();
        double *ms_3d_new_ptr = new double[size];
      
        for(int i = 0; i < size; i++)
            ms_3d_new_ptr[i] = (double) ms_3d_new.v()[i];
       
        // Recovering the dimensions of the images
        int width = hyperDataUpsampled.w();
        int height = hyperDataUpsampled.h();
        int hs_channels = hyperDataUpsampled.c();
        
        int width_new = im_3d_new.w();
        int height_new = im_3d_new.h();
        int ms_channels = ms_image.c();
        
        // Adapt filter parameters to size of comparison window
        int patch_size_2d = (2 * spat_patch_step + 1) * (2*spat_patch_step+1); // the size of the 2d patch
        double hSim2_2d = hSim * hSim * ((double) patch_size_2d);
        double hClose2 = hClose * hClose;
        
        std::vector<double>  diff_2d(patch_size_2d);
        std::vector<double>  diff_2_2d(patch_size_2d);
        std::vector<double>  tmp_pix_central_2d(patch_size_2d);
        std::vector<double>  tmp_pix_ser_2d(patch_size_2d);
        
        for(int h = 0; h < hs_channels; h++)
        {
            int lh = h * width * height;
            
            // Browsing the rows
            // The height here is the NUMBER OF ROWS
            for(int x = spat_patch_step; x < height_new - spat_patch_step; x++)
            {
                int lx = (x - spat_patch_step) * width;
                
                // Restraining the search in the spatial (x) domain in case we are near the borders
                int bound_start_x = MAX(x - reswind, spat_patch_step);
                int bound_end_x=MIN(x + reswind, height_new - spat_patch_step - 1);
                
                // Browsing the columns
                // The width here is THE NUMBER OF COLUMNS
                for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
                {
                    // Linear index of the current pixel at the center of the search window
                    int l_comp = (y - spat_patch_step) + lx; // Index for the comparison matrix
                    
                    // Restraining the search in the spatial (x) domain in case we are near the borders
                    int bound_start_y = MAX(y - reswind, spat_patch_step);
                    int bound_end_y = MIN(y + reswind, width_new - spat_patch_step - 1);
                    
                    // Declaring the vector where 3D unclassified weights computed with all
                    // The spectral weighted bands of the multispectral are going to be stored
                    std::vector<double> vect_xy_ms((2*reswind+1)*(2*reswind+1));
                    
                    // Computing the unclassified 3D weights on the multispectral image
                    nlweights_3Dall_ms_Smatrix_unclassified_double(ms_3d_new_ptr,vect_xy_ms, hSim, hClose, reswind, spat_patch_step,
                                                                   ms_channels, width_new, height_new, Smatrix, x, y, h);
                    
                    // Slicing (isolating) the 2D patch (centered on the central pixel (x,y,0)) out of the image
                    array_slicing_double(x-spat_patch_step, y-spat_patch_step, h, x+spat_patch_step, y+spat_patch_step,
                                         h, height_new, width_new, im_3d_new_ptr, tmp_pix_central_2d);
                    
                    // (Re)initializing variables
                    int w = 0; // browses the search window for each central pixel
                    int wcentral = 0;
                    double dMin_2d = dLarge;
                    double dSim_2d = 0e0;
                    double dClose = 0e0;
                    
                    // Auxiliary vector
                    std::vector<double> tmp_wxy((2*reswind+1)*(2*reswind+1));
                    std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));
                    
                    // Entering the search window
                    for(int x_ser = bound_start_x; x_ser <= bound_end_x; x_ser++) // x_ser (x in search window)
                    {
                        int lx_ser = (x_ser - spat_patch_step) * width;
                        
                        for(int y_ser = bound_start_y; y_ser <= bound_end_y; y_ser++) // y_ser (y in search window)
                        {
                            // Linear index of the current pixel in the search window
                            int l_ser_2d = (y_ser - spat_patch_step) + lx_ser;
                            
                            // Slicing (isolating) the 2D patch out of the image
                            array_slicing_double(x_ser-spat_patch_step, y_ser-spat_patch_step, h, x_ser + spat_patch_step,
                                                 y_ser + spat_patch_step, h, height_new, width_new, im_3d_new_ptr, tmp_pix_ser_2d);
                            
                            // BILATERAL DISTANCE WITH THE 2D PATCHS
                            // Computing the patch-based similarity distance
                            substract_two_vectors_double(tmp_pix_central_2d, tmp_pix_ser_2d, patch_size_2d, diff_2d);
                            power_two_vector_double(diff_2d, patch_size_2d, diff_2_2d);
                            dSim_2d=sum_values_vector_double(diff_2_2d, patch_size_2d);
                            dSim_2d /= hSim2_2d;
                            
                            // Computing the spatial distance
                            dClose = (double) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
                            dClose /= hClose2;
                            
                            // Storing the bilateral distances computed with the 3d patches and 2d patches
                            double d_bilateral_2d = dSim_2d + dClose;
                            tmp_wxy[w] = d_bilateral_2d;
                            
                            // Save w-position of central pixel in neighborhood and the minimum bilateral distance
                            if((x==x_ser) && (y==y_ser))
                            {
                                wcentral = w;
                            }
                            else
                            {
                                if(d_bilateral_2d < dMin_2d)
                                    dMin_2d = d_bilateral_2d;
                            }
                            
                            // Save neighbouring pixel position
                            tmp_pxy[w] = l_ser_2d;
                            
                            // Update index
                            w++;
                        }
                    }
                    
                    // Restraining the size of tmp_wxy to w
                    tmp_wxy.resize(w);
                    tmp_pxy.resize(w);
                    
                    // Assign minimum distance to central pixel
                    tmp_wxy[wcentral]=dMin_2d;
                    
                    // Adapt number of neighboring pixels to window size
                    int numSimPixels0 = MIN(numSimPixels, w);
                    sort_wrt_comp_vect_double(tmp_wxy, tmp_pxy, vect_xy_ms);
                    
                    // Keep only the first numSimPixels0 from wxy and posxy
                    tmp_wxy.resize(numSimPixels0);
                    tmp_pxy.resize(numSimPixels0);
                    
                    // Applying the exponential
                    double sum_weights = 0e0;
                    
                    for(int r = 0; r < numSimPixels0; r++)
                    {
                        double weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
                        tmp_wxy[r] = weight;
                        sum_weights += weight;
                    }
                    
                    // Normalization and applying the square root (sqrt)
                    if(flagNormalize)
                    {
                        if(sum_weights > dTiny)
                        {
                            for(int r = 0; r < numSimPixels0; r++)
                            {
                                double weight = (double) (tmp_wxy[r] / sum_weights);
                                tmp_wxy[r] = sqrt(weight);
                            }
                            
                        } else
                        {
                            tmp_wxy.clear();
                            tmp_pxy.clear();
                            
                            tmp_wxy.push_back(1e0);
                            tmp_pxy.push_back(l_comp);
                        }
                        
                    } else
                    {
                        for(int r = 0; r < numSimPixels0; r++)
                        {
                            double weight = (double) tmp_wxy[r];
                            tmp_wxy[r] = sqrt(weight);
                        }
                    }
                    
                    wxy.push_back(tmp_wxy);
                    pxy.push_back(tmp_pxy);
                }
            }
            
            // Looking for pyx, wyx and posw
            for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
            {
                int lx = (x - spat_patch_step) * width;
                
                // Restraining the search in the spatial (x) domain in case we are near the borders
                int bound_start_x = MAX(x - reswind, spat_patch_step);
                int bound_end_x = MIN(x + reswind, height_new - spat_patch_step - 1);
                
                // Browsing the columns
                // The width here is THE NUMBER OF COLUMNS
                for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
                {
                    // Linear index of the current pixel at the center of the search window
                    int l_2d = (y - spat_patch_step) + lx;
                    
                    // Restraining the search in the spatial (x) domain in case we are near the borders
                    int bound_start_y = MAX(y - reswind, spat_patch_step);
                    int bound_end_y = MIN(y + reswind, width_new - spat_patch_step - 1);
                    
                    // Auxiliary vector
                    std::vector<double> tmp_wyx((2*reswind+1)*(2*reswind+1));
                    std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
                    std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));
                    
                    // (Re)initializing variables
                    int w = 0;
                    
                    // Entering the search window
                    for(int x_ser = bound_start_x; x_ser <= bound_end_x; x_ser++) // x_ser (x in search window)
                    {
                        int lx_ser = (x_ser - spat_patch_step) * width;
                        
                        for(int y_ser = bound_start_y; y_ser <= bound_end_y; y_ser++) // y_ser (y in search window)
                        {
                            // The linear index of the current pixel in the search window
                            int l_ser = (y_ser - spat_patch_step) + lx_ser + lh;
                            int l_ser_2d = (y_ser - spat_patch_step) + lx_ser;
                            
                            // We verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
                            int pos = std::find(pxy[l_ser].begin(), pxy[l_ser].end(),l_2d)-pxy[l_ser].begin();
                            
                            // In the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
                            if(pos < (int) pxy[l_ser].size())
                            {
                                tmp_wyx[w] = wxy[l_ser][pos];
                                tmp_pyx[w] = l_ser_2d;
                                tmp_posw[w] = pos;
                                w++;
                            }
                        }
                    }
                    tmp_wyx.resize(w);
                    tmp_pyx.resize(w);
                    tmp_posw.resize(w);
                    
                    wyx.push_back(tmp_wyx);
                    pyx.push_back(tmp_pyx);
                    posw.push_back(tmp_posw);
                }
            }
        }
    }


    
    // Compute 3D distances on the multispectral image weighted with coefficient of the Smatrix.
    // The distances are unclassified.
    // Width_ms_new and height_ms_new are the spatial dimensions of the multispectral image after spatial extension.
    void nlweights_3Dall_ms_Smatrix_unclassified_double(double *ms_image_new_ptr, std::vector<double> &vect_xy, double hSim, double hClose,
                                                        int reswind, int spat_patch_step, int ms_channels, int width_ms_new, int height_ms_new,
                                                        double **Smatrix, int x, int y, int num_hs_channel)
    {
        // Adapt filter parameters to size of comparison window
        int patch_size_3d = (2 * spat_patch_step + 1) * (2 * spat_patch_step + 1) * (ms_channels); // the size of the 3D patch
        int dim = (2 * spat_patch_step + 1) * (2 * spat_patch_step + 1);
        double hSim2_3d = hSim * hSim * ((double) patch_size_3d);
        double hClose2 = hClose * hClose;
        
        std::vector<double> diff_3d(patch_size_3d);
        std::vector<double> diff_2_3d(patch_size_3d);
        std::vector<double> tmp_pix_central_3d(patch_size_3d);
        
        // Restraining the search in the spatial (x) domain in case we are near the borders
        int bound_start_x = MAX(x - reswind, spat_patch_step);
        int bound_end_x = MIN(x + reswind, height_ms_new - spat_patch_step - 1);
        
        // Restraining the search in the spatial (y) domain in case we are near the borders
        int bound_start_y = MAX(y - reswind, spat_patch_step);
        int bound_end_y = MIN(y + reswind, width_ms_new - spat_patch_step - 1);
        
        // Slicing (isolating) the 3D patch out of the image
        double *coef_tab = new double[ms_channels];
        double sum_coef = 0e0;
        
        for(int m = 0; m < ms_channels; m++)
        {
            double sval = Smatrix[m][num_hs_channel];
            coef_tab[m] = sval;
            sum_coef = sum_coef + sval;
        }
        
        int ind_cent = 0;
        
        for(int m = 0; m < ms_channels; m++)
        {
            double coef = coef_tab[m];
            array_slicing_coef_double(x - spat_patch_step, y - spat_patch_step, m, x + spat_patch_step, y + spat_patch_step, m,
                                      height_ms_new, width_ms_new, ms_image_new_ptr, ind_cent, tmp_pix_central_3d, coef);
            ind_cent = ind_cent + dim;
        }
        
        // (Re)initializing variables
        int w = 0; // browses the search window for each central pixel
        int wcentral = 0;
        double dMin_3d = dLarge;
        double dSim_3d = 0e0;
        double dClose = 0e0;
        
        // Entering the search window
        for(int x_ser = bound_start_x; x_ser <= bound_end_x; x_ser++) // x_ser (x in search window)
        {
            for(int y_ser = bound_start_y; y_ser <= bound_end_y; y_ser++) // y_ser (y in search window)
            {
                std::vector<double> tmp_pix_ser_3d(patch_size_3d);
                
                // Slicing (isolating) the 3D patch out of the image
                int ind_ser = 0;
                
                for(int m = 0; m < ms_channels; m++)
                {
                    double coef = coef_tab[m];
                    array_slicing_coef_double(x_ser-spat_patch_step, y_ser-spat_patch_step, m, x_ser+spat_patch_step, y_ser+spat_patch_step, m,
                                              height_ms_new, width_ms_new, ms_image_new_ptr, ind_ser, tmp_pix_ser_3d, coef);
                    ind_ser = ind_ser + dim;
                }
                
                // BILATERAL DISTANCE WITH THE 3D PATCHS
                // Computing the patch-based similarity distance
                substract_two_vectors_double(tmp_pix_central_3d, tmp_pix_ser_3d, patch_size_3d, diff_3d);
                power_two_vector_double(diff_3d, patch_size_3d, diff_2_3d);
                dSim_3d = sum_values_vector_double(diff_2_3d, patch_size_3d);
                
                double hSim2 = hSim2_3d * sum_coef;
                dSim_3d /= hSim2;
                
                // Computing the spatial distance
                dClose=(double) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
                dClose /= hClose2;
                
                // Storing the bilateral distances computed with the 3D patches and 2D patches
                double d_bilateral_3d = dSim_3d + dClose;
                
                vect_xy[w] = d_bilateral_3d;
                
                // Save w-position of central pixel in neighborhood and the minimum bilateral distance
                if((x == x_ser) && (y == y_ser))
                {
                    wcentral = w;
                }
                else
                {
                    if(d_bilateral_3d < dMin_3d)
                        dMin_3d = d_bilateral_3d;
                }
                
                // Update index
                w++;
            }
        }
        
        // Restraining the size of WandPxy[l] to w
        vect_xy.resize(w);
        
        // Assign minimum distance to central pixel
        vect_xy[wcentral] = dMin_3d;
    }



    // Sort with respect to a comparison vector (ascending order)
    void sort_wrt_comp_vect_double(std::vector<double> &values_vect, std::vector<int> &positions_vect, const std::vector<double> &comparison_vect)
    {
        int size = (int) values_vect.size();
        std::vector<ValPosCompVect> vector(size);
        
        // Copy everything in a structure
        for(int i = 0; i < size; i++)
        {
            vector[i].value = values_vect[i];
            vector[i].comp_value = comparison_vect[i];
            vector[i].position = positions_vect[i];
        }
        
        // Sorting (in the ascending direction)
        std::sort(vector.begin(), vector.end(), [](const ValPosCompVect &a, const ValPosCompVect &b) { return a.comp_value < b.comp_value; } );
        
        // Recopy everything into the original vectors
        for(int i = 0; i < size; i++)
        {
            values_vect[i] = vector[i].value;
            //comparison_vect[i] = vector[i].comp_value;
            positions_vect[i] = vector[i].position;
        }
    }
    
    
    
    void array_slicing_coef_double(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height,
                                   int width, double *im_to_slice, int ind, std::vector<double> &sliced_vector, double coef)
    {
        int cpt = ind;
        
        for(int h = h_start; h <= h_end; h++)
        {
            for(int i = i_start; i <= i_end; i++)
            {
                for(int j = j_start; j <= j_end; j++)
                {
                    sliced_vector[cpt] = coef * im_to_slice[j+(i*width)+(h*height*width)];
                    cpt++;
                }
            }
        }
    }

    
    
    void array_slicing_double(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height,
                              int width, double *im_to_slice, std::vector<double> &sliced_vector)
    {
        int cpt = 0;
        
        for(int h = h_start; h <= h_end; h++)
        {
            for(int i = i_start; i <= i_end; i++)
            {
                for(int j = j_start; j <= j_end; j++)
                {
                    sliced_vector[cpt] = im_to_slice[j+(i*width)+(h*height*width)];
                    cpt++;
                }
            }
        }
    }

    
    
    double sum_values_vector_double(std::vector<double> &vec, int size_vec)
    {
        double sum = 0e0;
        
        for(int j = 0; j < size_vec; j++)
            sum = sum + vec[j];
        
        return sum;
    }



    void substract_two_vectors_double(std::vector<double> &vec1, std::vector<double> &vec2, int size_vec, std::vector<double> &res_vector)
    {
        
        for(int i = 0; i < size_vec; i++)
            res_vector[i] = vec1[i] - vec2[i];
    }


    
    void power_two_vector_double(std::vector<double> &vector, int size_vec, std::vector<double> &res_vector)
    {
        for(int k = 0; k < size_vec; k++)
            res_vector[k] = vector[k] * vector[k];
    }


// Extends the SPATIAL dimensions of an image by mirror effect which helps in case the treatment of the image takes the borders into
// account. The image to extend should be a 3D image.
libUSTG::cflimage extend_image_dimensions_spat(libUSTG::cflimage im_to_extend, int spat_patch_step)
{
	/* the image to extend is the image from the cflimage class
	 * so that we can access the dimensions of the image easily
	 */

	/*extracting the dimensions of the image to extend*/

	int nr_old=im_to_extend.h(); //number of rows of the image to extend
	int nc_old=im_to_extend.w(); //number of columns of the image to extend
	int b=im_to_extend.c(); //number of spectral bands of the image to extend (it's not going to change)

	/* peut etre qu'il faut faire un test entre le pas du patch et le nombre de lignes
	 * et de colonnes de l'image à etendre, si le pas du patch est supérieur au nombre de lignes et
	 * de colonnes de l'image à étendre je trouve que cela n'a pas de sens, si ??*/

	/*Test in case the patch step is bigger than the number of rows and cols of the image to extend
	 *If it's the case (which I don't find interesting) and the image is not extended.
	 */

	if((spat_patch_step>nr_old)&&(spat_patch_step>nc_old))
	{
		printf("Error, the patch step is bigger than the rows and the columns if the image to extend, not an interesting case \n");
		return im_to_extend;
	}


	/*Test in case the patch step is null*/
	if(spat_patch_step==0)
	{
		return im_to_extend;
	}


	/* step1: creating a matrix with the new spatial dimensions */
	int nr_new=nr_old+(2*spat_patch_step);
	int nc_new=nc_old+(2*spat_patch_step);

	libUSTG::cflimage extended_image;
	extended_image.create(nc_new,nr_new,b);

	int k,j,i,m=0;
	for(k=0;k<b;k++)
	{

		/*step2: copy of the content of the old matrix into the new matrix*/
		/*the browsing is done in rows*/
		for(j=spat_patch_step;j<=(nc_new-spat_patch_step-1);j++)
		{
			for(i=spat_patch_step;i<=(nr_new-spat_patch_step-1);i++)
			{
				extended_image[(k*nc_new*nr_new)+(j*nr_new)+i]=im_to_extend[m];
				m++;

			}
		}

		/*step3: filling the columns*/
		for(i=spat_patch_step;i<=(nr_new-spat_patch_step-1);i++)
		{
			/*left columns*/
			int j_prime=(2*spat_patch_step)-1;
			for(j=0;j<=(spat_patch_step-1);j++)
			{
				extended_image[(k*nc_new*nr_new)+(j*nr_new)+i]=extended_image[(k*nc_new*nr_new)+(j_prime*nr_new)+i];
				j_prime--;
			}

			/*right columns*/
			int j_prime_second=nc_new-(2*spat_patch_step);
			for(j=nc_new-1;j>=(nc_new-spat_patch_step);j--)
			{
				extended_image[(k*nc_new*nr_new)+(j*nr_new)+i]=extended_image[(k*nc_new*nr_new)+(j_prime_second*nr_new)+i];
				j_prime_second++;
			}
		}

		/*step4: filling the rows*/
		for(j=0;j<nc_new;j++)
		{
			/*above rows*/
			int i_prime=(2*spat_patch_step)-1;
			for(i=0;i<=(spat_patch_step-1);i++)
			{
				extended_image[(k*nc_new*nr_new)+(j*nr_new)+i]=extended_image[(k*nc_new*nr_new)+(j*nr_new)+i_prime];
				i_prime--;
			}

			/*below rows*/
			int i_prime_second=nr_new-(2*spat_patch_step);
			for(i=nr_new-1;i>=(nr_new-spat_patch_step);i--)
			{
				extended_image[(k*nc_new*nr_new)+(j*nr_new)+i]=extended_image[(k*nc_new*nr_new)+(j*nr_new)+i_prime_second];
				i_prime_second++;
			}
		}

	}


	return extended_image;

}


} // libUSTGHYPERDOUBLE
