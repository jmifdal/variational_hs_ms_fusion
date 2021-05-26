include "libHyperspectralWeights.h"


namespace libUSTGHYPER
{
// Computes the sqrt(weight) on the MS data, so 3D patches are considered
// The patch distance across channels is coupled with the L2 norm.
// The same weight distribution is obtained for all HS channels.
// Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
void nlweights_ms(float **mData, std::vector< std::vector<float> > &wxy, std::vector< std::vector<int> > &posxy,
		std::vector< std::vector<float> > &wyx, std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw,
		float hPatch, float hSpatial, int numNeighbours, int reswind, int compwind, int flagNormalize, int ms_channels,
		int width, int height)
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

	// Filtering parameter for spatial distance
	float hSpatial2 = hSpatial * hSpatial;

	// For each pixel l, compute weights w(l,l0)
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
			float hPatch2 = hPatch * hPatch * compdim;

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
					float dSim = libUSTG::fiL2FloatDist(mData, mData, x, y, i, j, compwind0, compwind0, ms_channels, width, width);
					dSim /= hPatch2;

					// Compute spatial distance
					float dClose = (float) ((x - i) * (x - i) + (y - j) * (y - j));
					dClose /= hSpatial2;

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


	// For each pixel l, compute weights w(l0,l)
	for(int y = 0; y < height; y++)
	{
		// Fix row
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

			wyx.push_back(auxw);
			posyx.push_back(auxp);
			posw.push_back(auxpw);
		}
	}

	// Delete allocated memory
	delete[] lut;
	delete[] dist;
	delete[] pos;
}

// Computes the sqrt(weight) on the MS data, so 3D patches are considered.
    // The contribution of each MS channel to the patch-based distance is weighted by the coefficients of the spectral downsampling matrix.
    // A different weight distribution is obtained for each HS channel.
    // Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
void nlweights_ms_Smatrix_new(libUSTG::cflimage ms_image,float **Smatrix ,
		std::vector< std::vector<float> > &wxy,std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
		std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, int hs_channels, float hSim,
		float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize)
{

	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	//extension of the images
	//in this case the image should be extended only spatially
	libUSTG::cflimage ms_3d_new=extend_image_dimensions_spat(ms_image,spat_patch_step); //spatial extension
	float *ms_3d_new_ptr=ms_3d_new.d_v;

	//recovering the dimensions of the images
	int width_ms=ms_image.w();
	int height_ms=ms_image.h();

	int width_ms_new=ms_3d_new.w();
	int height_ms_new=ms_3d_new.h();

	int ms_channels=ms_image.c();

	// Adapt filter parameters to size of comparison window
	int patch_size_3d=(2 * spat_patch_step + 1)*(2*spat_patch_step+1)*ms_channels; //The size of the 3d patch
	int dim =(2 * spat_patch_step + 1)*(2*spat_patch_step+1);//the size of the 2d patch

	float hSim2_3d = hSim * hSim * patch_size_3d;

	float hClose2 = hClose * hClose;

	std::vector<float>  diff_3d(patch_size_3d);
	std::vector<float>  diff_2_3d(patch_size_3d);

	std::vector<float>  tmp_pix_central_3d(patch_size_3d);
	std::vector<float>  tmp_pix_ser_3d(patch_size_3d);

	//slicing (isolating) the 3d patch out of the image
	for(int h=0;h<hs_channels;h++)
	{

		int lh = h*width_ms*height_ms;

		float *coef_tab=new float[ms_channels];
		float sum_coef=0.0;
		for(int m=0;m<ms_channels;m++)
		{
			float sval=Smatrix[m][h];
			coef_tab[m]=sval;
			sum_coef=sum_coef+sval;
		}


		//browsing the rows
		//the height here is the NUMBER OF ROWS
		for(int x = spat_patch_step; x < height_ms_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width_ms;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_ms_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_ms_new-spat_patch_step; y++)
			{
				//linear index of the current pixel at the center of the search window
				int l_comp=(y-spat_patch_step)+lx;//index for the comparison matrix

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_ms_new-spat_patch_step-1);

				int ind_cent=0;
				for(int m=0;m<ms_channels;m++)
				{
					float coef=coef_tab[m];
					array_slicing_coef(x-spat_patch_step,y-spat_patch_step,m,x+spat_patch_step
							,y+spat_patch_step,m,width_ms_new,height_ms_new,ms_3d_new_ptr,ind_cent,tmp_pix_central_3d,coef);
					ind_cent=ind_cent+dim;
				}

				//(re)initializing variables
				int w=0; //browses the search window for each central pixel
				int wcentral=0;
				float dMin_3d = fLarge;
				float dSim_3d=0;
				float dClose=0;

				//auxiliary vector
				std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width_ms;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//linear index of the current pixel in the search window
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						//slicing (isolating) the 3d patch out of the image
						int ind_ser=0;
						for(int m=0;m<ms_channels;m++)
						{
							float coef=coef_tab[m];
							array_slicing_coef(x_ser-spat_patch_step,y_ser-spat_patch_step,m,x_ser+spat_patch_step
									,y_ser+spat_patch_step,m,height_ms_new,width_ms_new,ms_3d_new_ptr,ind_ser,tmp_pix_ser_3d,coef);
							ind_ser=ind_ser+dim;
						}

						//BILATERAL DISTANCE WITH THE 3D PATCHS
						//computing the patch-based similarity distance
						substract_two_vectors(tmp_pix_central_3d,tmp_pix_ser_3d,patch_size_3d,diff_3d);

						power_two_vector(diff_3d,patch_size_3d,diff_2_3d);

						dSim_3d=sum_values_vector(diff_2_3d,patch_size_3d);

						float hSim2=hSim2_3d*sum_coef;

						dSim_3d /= hSim2;

						//computing the spatial distance
						dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
						dClose /= hClose2;

						//storing the bilateral distances computed with the 3d patches and 2d patches
						float d_bilateral_3d=dSim_3d+dClose;
						tmp_wxy[w]=d_bilateral_3d;

						// Save w-position of central pixel in neighborhood and the minimum bilateral distance
						if((x==x_ser) && (y==y_ser))
						{
							wcentral = w;
						}
						else
						{
							if(d_bilateral_3d < dMin_3d)
								dMin_3d = d_bilateral_3d;
						}

						//Save neighbouring pixel position
						tmp_pxy[w]=l_ser_2d;

						// Update index
						w++;
					}
				}

				//restraining the size of tmp_wxy to w
				tmp_wxy.resize(w);
				tmp_pxy.resize(w);

				//Assign minimum distance to central pixel
				tmp_wxy[wcentral]=dMin_3d;

				//Adapt number of neighboring pixels to window size
				int numSimPixels0 = MIN(numSimPixels, w);

				//sorting the weights and their positions in an ascending order
				sort_vects(tmp_wxy, tmp_pxy);

				//Keep only the first numSimPixels0 from wxy and posxy
				tmp_wxy.resize(numSimPixels0);
				tmp_pxy.resize(numSimPixels0);

				// applying the exponential
				float sum_weights=0.0f;
				for(int r = 0; r < numSimPixels0; r++)
				{
					float weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
					tmp_wxy[r]=weight;
					sum_weights += weight;
				}

				//normalization and applying the square root (sqrt)
				if(flagNormalize)
				{
					if(sum_weights > fTiny)
					{
						for(int r = 0; r < numSimPixels0; r++)
						{
							float weight = (float) tmp_wxy[r] / sum_weights;
							tmp_wxy[r] = sqrt(weight);
							//std::cout<<"weight = "<<tmp_wxy[r]<<std::endl;
						}

					} else
					{
						tmp_wxy.clear();
						tmp_pxy.clear();

						tmp_wxy.push_back(1.0f);
						tmp_pxy.push_back(l_comp);
					}

				} else
				{
					for(int r = 0; r < numSimPixels0; r++)
					{
						float weight = (float) tmp_wxy[r];
						tmp_wxy[r] = sqrt(weight);
					}
				}

				wxy.push_back(tmp_wxy);
				pxy.push_back(tmp_pxy);
			}
		}

		//looking for pyx, wyx and posw
		for(int x = spat_patch_step; x < height_ms_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width_ms;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_ms_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_ms_new-spat_patch_step; y++)
			{

				//linear index of the current pixel at the center of the search window
				int l_2d=(y-spat_patch_step)+lx;

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_ms_new-spat_patch_step-1);

				//auxiliary vector
				std::vector<float> tmp_wyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));

				//(re)initializing variables
				int w=0;

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width_ms;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//the linear index of the current pixel in the search window
						int l_ser=(y_ser-spat_patch_step)+lx_ser+lh;
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						//we verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
						int pos=std::find(pxy[l_ser].begin(), pxy[l_ser].end(),l_2d)-pxy[l_ser].begin();

						//in the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
						if( pos < (int) pxy[l_ser].size() )
						{
							tmp_wyx[w]=wxy[l_ser][pos];
							tmp_pyx[w]=l_ser_2d;
							tmp_posw[w]=pos;
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


// Computes the sqrt(weight) on the MS data, so 3D patches are considered.
// The contribution of each MS channel to the patch-based distance is weighted by the coefficients of the spectral downsampling matrix.
// A different weight distribution is obtained for each HS channel.
// Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
void nlweights_ms_Smatrix(float **mData, float **Smatrix, std::vector< std::vector< std::vector<float> > > &wxy,
		std::vector< std::vector< std::vector<int> > > &posxy, std::vector< std::vector< std::vector<float> > > &wyx,
		std::vector< std::vector< std::vector<int> > > &posyx, std::vector< std::vector< std::vector<int> > > &posw,
		float hPatch, float hSpatial, int numNeighbours, int reswind, int compwind, int flagNormalize,
		int hs_channels, int ms_channels, int width, int height)
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

	// Filtering parameter for spatial distance
	float hSpatial2 = hSpatial * hSpatial;

	// For each hyperspectral channel
	for(int h = 0; h < hs_channels; h++)
	{
		{
			// Auxiliar vectors
			std::vector< std::vector<float> > wxyaux;
			std::vector< std::vector<int> > posxyaux;

			// For each pixel l, compute weights w(l,l0)
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
					float hPatch2 = hPatch * hPatch * compdim;

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
							float dSim = 0.0f;
							float sumCoef = 0.0f;

							for(int m = 0; m < ms_channels; m++)
							{
								float coef = Smatrix[m][h];
								dSim += (coef * libUSTG::fiL2FloatDist(mData[m], mData[m], x, y, i, j, compwind0, compwind0, width, width));
								sumCoef += coef;
							}

							dSim /= (hPatch2 * sumCoef * (float) ms_channels);

							// Compute spatial distance
							float dClose = (float) ((x - i) * (x - i) + (y - j) * (y - j));
							dClose /= hSpatial2;

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

					wxyaux.push_back(auxw);
					posxyaux.push_back(auxp);
				}
			}

			//Save outputs
			wxy.push_back(wxyaux);
			posxy.push_back(posxyaux);
		}

		{
			// Auxiliar vectors
			std::vector< std::vector<float> > wyxaux;
			std::vector< std::vector<int> > posyxaux;
			std::vector< std::vector<int> > poswaux;

			// For each pixel l, compute weights w(l0,l)
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
							for(int w = 0; w < (int) wxy[h][l0].size(); w++)
							{
								if(posxy[h][l0][w] == l)
								{
									float weight = wxy[h][l0][w];

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

					wyxaux.push_back(auxw);
					posyxaux.push_back(auxp);
					poswaux.push_back(auxpw);
				}
			}

			// Save outputs
			wyx.push_back(wyxaux);
			posyx.push_back(posyxaux);
			posw.push_back(poswaux);
		}
	}

	// Delete allocated memory
	delete[] lut;
	delete[] dist;
	delete[] pos;
}



// For each HS channel, computes the sqrt(weight) on the corresponding upsampled HS image, so 2D patches are considered.
// A different weight distribution is obtained for each HS channel.
// Weight structure: wxy[h][i][n] for each channel h, each pixel i, and each neighbour n.
void nlweights_hs_band(float **hDataInt, std::vector< std::vector< std::vector<float> > > &wxy,
		std::vector< std::vector< std::vector<int> > > &posxy, std::vector< std::vector< std::vector<float> > > &wyx,
		std::vector< std::vector< std::vector<int> > > &posyx, std::vector< std::vector< std::vector<int> > > &posw,
		float hPatch, float hSpatial, int numNeighbours, int reswind, int compwind, int flagNormalize, int hs_channels,
		int width, int height)
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

	// Filtering parameter for spatial distance
	float hSpatial2 = hSpatial * hSpatial;

	// For each hyperspectral channel
	for(int h = 0; h < hs_channels; h++)
	{
		{
			// Auxiliar vectors
			std::vector< std::vector<float> > wxyaux;
			std::vector< std::vector<int> > posxyaux;

			// For each pixel l, compute weights w(l,l0)
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
					float hPatch2 = hPatch * hPatch * compdim;

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
							float dSim = libUSTG::fiL2FloatDist(hDataInt[h], hDataInt[h], x, y, i, j, compwind0, compwind0, width, width);
							dSim /= hPatch2;

							// Compute spatial distance
							float dClose = (float) ((x - i) * (x - i) + (y - j) * (y - j));
							dClose /= hSpatial2;

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

					wxyaux.push_back(auxw);
					posxyaux.push_back(auxp);
				}
			}

			//Save outputs
			wxy.push_back(wxyaux);
			posxy.push_back(posxyaux);
		}

		{
			// Auxiliar vectors
			std::vector< std::vector<float> > wyxaux;
			std::vector< std::vector<int> > posyxaux;
			std::vector< std::vector<int> > poswaux;

			// For each pixel l, compute weights w(l0,l)
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
							for(int w = 0; w < (int) wxy[h][l0].size(); w++)
							{
								if(posxy[h][l0][w] == l)
								{
									float weight = wxy[h][l0][w];

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

					wyxaux.push_back(auxw);
					posyxaux.push_back(auxp);
					poswaux.push_back(auxpw);
				}
			}

			// Save outputs
			wyx.push_back(wyxaux);
			posyx.push_back(posyxaux);
			posw.push_back(poswaux);
		}
	}

	// Delete allocated memory
	delete[] lut;
	delete[] dist;
	delete[] pos;
}



// Computes the sqrt(weights) by first choosing most similar neighbours on the full MS image (3D patches) and computing then
// the weights using 2D patches centered at the pre-selected pixels on the corresponding band.
// A different weight distribution is obtained for each HS channel.
// Weight structure: wxy[h*dim+i][n] for each channel h, each pixel i, and each neighbour n.
void nlweights_ms_3Dclassify_2Dband(libUSTG::cflimage hyperDataUpsampled, libUSTG::cflimage ms_image,std::vector< std::vector<float> > &wxy,
		std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
		std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, float hSim,
		float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize)
{


	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	//extension of the image
	//in this case the image should be extended only spatially
	libUSTG::cflimage im_3d_new=extend_image_dimensions_spat(hyperDataUpsampled,spat_patch_step); //spatial extension
	float *im_3d_new_ptr=im_3d_new.d_v;

	libUSTG::cflimage ms_image_new=extend_image_dimensions_spat(ms_image,spat_patch_step); //spatial extension
	float *ms_image_new_ptr=ms_image_new.d_v;

	//recovering the dimensions of the images
	int width=hyperDataUpsampled.w();
	int height=hyperDataUpsampled.h();
	int hs_channels=hyperDataUpsampled.c();
	int ms_channels=ms_image.c();

	int width_new=im_3d_new.w();
	int height_new=im_3d_new.h();

	// Adapt filter parameters to size of comparison window
	int patch_size_2d =(2 * spat_patch_step + 1)*(2*spat_patch_step+1);//the size of the 2d patch

	float hSim2_2d = hSim * hSim * patch_size_2d;

	float hClose2 = hClose * hClose;

	std::vector<float>  diff_2d(patch_size_2d);
	std::vector<float>  diff_2_2d(patch_size_2d);

	std::vector<float>  tmp_pix_central_2d(patch_size_2d);
	std::vector<float>  tmp_pix_ser_2d(patch_size_2d);

	//declaring the matrix where 3d unclassified weights computed with all the spectral bands (of the multispectral)
	//are going to be stored
	std::vector< std::vector<float> > wxy_3D_ms;

	//Computing the unclassified 3d weights on the multispectral image
	nlweights_3Dall_unclassified(ms_image_new_ptr, wxy_3D_ms, hSim, hClose, reswind, spat_patch_step, ms_channels, width_new,
			height_new);

	for(int h=0;h<hs_channels;h++)
	{
		//std::cout<<" channel "<<h<<std::endl;
		int lh = h*width*height;
		//browsing the rows
		//the height here is the NUMBER OF ROWS
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			{
				//linear index of the current pixel at the center of the search window
				int l_comp=(y-spat_patch_step)+lx;//index for the comparison matrix

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//slicing (isolating) the 2d patch (centered on the central pixel (x,y,0)) out of the image
				array_slicing(x-spat_patch_step,y-spat_patch_step,h,x+spat_patch_step
						,y+spat_patch_step,h,height_new,width_new,im_3d_new_ptr,tmp_pix_central_2d);

				//(re)initializing variables
				int w=0; //browses the search window for each central pixel
				int wcentral=0;
				float dMin_2d = fLarge;
				float dSim_2d=0;
				float dClose=0;

				//auxiliary vector
				std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//linear index of the current pixel in the search window
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;


						//slicing (isolating) the 2d patch out of the image
						array_slicing(x_ser-spat_patch_step,y_ser-spat_patch_step,h,x_ser+spat_patch_step,
								y_ser+spat_patch_step,h,height_new,width_new,im_3d_new_ptr,tmp_pix_ser_2d);

						//BILATERAL DISTANCE WITH THE 2D PATCHS
						//computing the patch-based similarity distance
						substract_two_vectors(tmp_pix_central_2d,tmp_pix_ser_2d,patch_size_2d,diff_2d);

						power_two_vector(diff_2d,patch_size_2d,diff_2_2d);

						dSim_2d=sum_values_vector(diff_2_2d,patch_size_2d);

						dSim_2d /= hSim2_2d;

						//computing the spatial distance
						dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
						dClose /= hClose2;


						//storing the bilateral distances computed with the 3d patches and 2d patches
						float d_bilateral_2d=dSim_2d+dClose;
						tmp_wxy[w]=d_bilateral_2d;

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

						//Save neighbouring pixel position
						tmp_pxy[w]=l_ser_2d;

						// Update index
						w++;


					}
				}

				//restraining the size of tmp_wxy to w
				tmp_wxy.resize(w);
				tmp_pxy.resize(w);


				//Assign minimum distance to central pixel
				tmp_wxy[wcentral]=dMin_2d;

				//Adapt number of neighboring pixels to window size
				int numSimPixels0 = MIN(numSimPixels, w);

				std::vector<float> vec=wxy_3D_ms[l_comp];
				sort_wrt_comp_vect(tmp_wxy, tmp_pxy, vec);

				//Keep only the first numSimPixels0 from wxy and posxy
				tmp_wxy.resize(numSimPixels0);
				tmp_pxy.resize(numSimPixels0);

				// applying the exponential
				float sum_weights=0.0f;
				for(int r = 0; r < numSimPixels0; r++)
				{
					float weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
					tmp_wxy[r]=weight;
					sum_weights += weight;
				}

				//normalization and applying the square root (sqrt)
				if(flagNormalize)
				{
					if(sum_weights > fTiny)
					{
						for(int r = 0; r < numSimPixels0; r++)
						{
							float weight = (float) tmp_wxy[r] / sum_weights;
							tmp_wxy[r] = sqrt(weight);
						}

					} else
					{
						tmp_wxy.clear();
						tmp_pxy.clear();

						tmp_wxy.push_back(1.0f);
						tmp_pxy.push_back(l_comp);
					}

				} else
				{
					for(int r = 0; r < numSimPixels0; r++)
					{
						float weight = (float) tmp_wxy[r];
						tmp_wxy[r] = sqrt(weight);
					}
				}

				wxy.push_back(tmp_wxy);
				pxy.push_back(tmp_pxy);
			}
		}

		//looking for pyx, wyx and posw
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			{

				//linear index of the current pixel at the center of the search window
				//int l=(y-spat_patch_step)+lx+lh;
				int l_2d=(y-spat_patch_step)+lx;

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//auxiliary vector
				std::vector<float> tmp_wyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));

				//(re)initializing variables
				int w=0;

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//the linear index of the current pixel in the search window
						int l_ser=(y_ser-spat_patch_step)+lx_ser+lh;
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						//we verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
						int pos=std::find(pxy[l_ser].begin(), pxy[l_ser].end(),l_2d)-pxy[l_ser].begin();

						//in the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
						if( pos < (int) pxy[l_ser].size() )
						{
							tmp_wyx[w]=wxy[l_ser][pos];
							tmp_pyx[w]=l_ser_2d;
							tmp_posw[w]=pos;
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





//In this function we classify 2D patch of the hyperspectral image based 3D
//patches computed on the multispectral image and weighted with the Smatrix

//TODO : maybe you should add tests to verify that the upsampled hyperspectral image and the multispectral image
//have the same spatial resolution. And also test if the number of the hyperspectral bands of the upsampled hyperspectral image and
//those of the multispectral bands are consistent of the dimensions of the Smatrix

void nlweights_ms_3Dclassify_Smatrix_2Dband(libUSTG::cflimage hyperDataUpsampled,libUSTG::cflimage ms_image,float **Smatrix ,
		std::vector< std::vector<float> > &wxy,std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
		std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, float hSim,
		float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize)
{

	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	//extension of the images
	//in this case the image should be extended only spatially
	libUSTG::cflimage im_3d_new=extend_image_dimensions_spat(hyperDataUpsampled,spat_patch_step); //spatial extension
	float *im_3d_new_ptr=im_3d_new.d_v;

	libUSTG::cflimage ms_3d_new=extend_image_dimensions_spat(ms_image,spat_patch_step); //spatial extension
	float *ms_3d_new_ptr=ms_3d_new.d_v;

	//recovering the dimensions of the images
	int width=hyperDataUpsampled.w();
	int height=hyperDataUpsampled.h();
	int hs_channels=hyperDataUpsampled.c();

	int width_new=im_3d_new.w();
	int height_new=im_3d_new.h();

	int ms_channels=ms_image.c();

	// Adapt filter parameters to size of comparison window
	int patch_size_2d =(2 * spat_patch_step + 1)*(2*spat_patch_step+1);//the size of the 2d patch

	float hSim2_2d = hSim * hSim * patch_size_2d;

	float hClose2 = hClose * hClose;

	std::vector<float>  diff_2d(patch_size_2d);
	std::vector<float>  diff_2_2d(patch_size_2d);

	std::vector<float>  tmp_pix_central_2d(patch_size_2d);
	std::vector<float>  tmp_pix_ser_2d(patch_size_2d);

	//DEBUGG
	//hs_channels=1;

	//DEBUGG
	//for(int h=64;h<65;h++)
	for(int h=0;h<hs_channels;h++)
	{
		int lh = h*width*height;

		//browsing the rows
		//the height here is the NUMBER OF ROWS
		//DEBUGG
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		//for(int x = 123+spat_patch_step; x <123+spat_patch_step+1 ; x++)
		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			//DEBUGG
			int cpt_norm=0;
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			//for(int y = 118+spat_patch_step; y < 118+spat_patch_step+1; y++)
			{
				//linear index of the current pixel at the center of the search window
				int l_comp=(y-spat_patch_step)+lx;//index for the comparison matrix


				//std::cout<<"h= "<<h<<", l_comp= "<<l_comp<<std::endl;		


				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);


				//declaring the vector where 3d unclassified weights computed with all
				//the spectral weighted bands of the multispectral are going to be stored
				std::vector<float> vect_xy_ms((2*reswind+1)*(2*reswind+1));

				//Computing the unclassified 3d weights on the multispectral image
				nlweights_3Dall_ms_Smatrix_unclassified(ms_3d_new_ptr,vect_xy_ms, hSim, hClose, reswind, spat_patch_step,
						ms_channels, width_new, height_new, Smatrix, x, y, h);

				/*for(int i=0; i<vect_xy_ms.size();i++)
				{
					std::cout<<"vect_xy_ms= "<<vect_xy_ms[i]<<std::endl;
				}*/


				//slicing (isolating) the 2d patch (centered on the central pixel (x,y,0)) out of the image
				array_slicing(x-spat_patch_step,y-spat_patch_step,h,x+spat_patch_step
						,y+spat_patch_step,h,height_new,width_new,im_3d_new_ptr,tmp_pix_central_2d);

				//(re)initializing variables
				int w=0; //browses the search window for each central pixel
				int wcentral=0;
				float dMin_2d = fLarge;
				float dSim_2d=0;
				float dClose=0;

				//auxiliary vector
				std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));

				//entering the search window

				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//linear index of the current pixel in the search window
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;


						//slicing (isolating) the 2d patch out of the image
						array_slicing(x_ser-spat_patch_step,y_ser-spat_patch_step,h,x_ser+spat_patch_step,
								y_ser+spat_patch_step,h,height_new,width_new,im_3d_new_ptr,tmp_pix_ser_2d);

						//BILATERAL DISTANCE WITH THE 2D PATCHS
						//computing the patch-based similarity distance
						substract_two_vectors(tmp_pix_central_2d,tmp_pix_ser_2d,patch_size_2d,diff_2d);

						power_two_vector(diff_2d,patch_size_2d,diff_2_2d);

						dSim_2d=sum_values_vector(diff_2_2d,patch_size_2d);

						dSim_2d /= hSim2_2d;

						//computing the spatial distance
						dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
						dClose /= hClose2;
						//printf("I'm here 3\n");

						//storing the bilateral distances computed with the 3d patches and 2d patches
						float d_bilateral_2d=dSim_2d+dClose;
						tmp_wxy[w]=d_bilateral_2d;

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

						//Save neighbouring pixel position
						tmp_pxy[w]=l_ser_2d;

						// Update index
						w++;


					}
				}

				//DEBUGG
				/*for(int i=0; i<tmp_wxy.size();i++)
                                {
                                        std::cout<<"tmp_wxy= "<<tmp_wxy[i]<<std::endl;
                                }*/


				//restraining the size of tmp_wxy to w
				tmp_wxy.resize(w);
				tmp_pxy.resize(w);

				//Assign minimum distance to central pixel
				tmp_wxy[wcentral]=dMin_2d;

				//Adapt number of neighboring pixels to window size
				int numSimPixels0 = MIN(numSimPixels, w);

				sort_wrt_comp_vect(tmp_wxy, tmp_pxy, vect_xy_ms);

				//Keep only the first numSimPixels0 from wxy and posxy
				tmp_wxy.resize(numSimPixels0);
				tmp_pxy.resize(numSimPixels0);


				// applying the exponential
				float sum_weights=0.0f;
				for(int r = 0; r < numSimPixels0; r++)
				{
					float weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
					tmp_wxy[r]=weight;
					sum_weights += weight;
				}


				/*for(int i=0; i<tmp_wxy.size();i++)
                                {
                                        std::cout<<"tmp_wxy= "<<tmp_wxy[i]<<", tmp_pxy= "<<tmp_pxy[i]<<std::endl;
                                }*/


				//normalization and applying the square root (sqrt)
				if(flagNormalize)
				{
				
					if(sum_weights > fTiny)
					{
						for(int r = 0; r < numSimPixels0; r++)
						{
							float weight = (float) tmp_wxy[r] / sum_weights;
							tmp_wxy[r] = sqrt(weight);
						}

					} else
					{
						printf("Am I here ? \n");
						tmp_wxy.clear();
						tmp_pxy.clear();

						tmp_wxy.push_back(1.0f);
						tmp_pxy.push_back(l_comp);

						//std::cout<<"l_comp= "<<l_comp<<", l_ser_2d= "<<l_ser_2d<<std::endl;
						//cpt_norm=cpt_norm+1;
					}

				} else
				{
					for(int r = 0; r < numSimPixels0; r++)
					{
						float weight = (float) tmp_wxy[r];
						tmp_wxy[r] = sqrt(weight);
					}
				}

				wxy.push_back(tmp_wxy);
				pxy.push_back(tmp_pxy);
				
				
				/*for(int i=0; i<tmp_wxy.size();i++)
                                {
                                        std::cout<<"(after norm)tmp_wxy= "<<tmp_wxy[i]<<", tmp_pxy= "<<tmp_pxy[i]<<std::endl;
                                }*/

			}
			//std::cout<<"cpt_norm= "<<cpt_norm<<std::endl;
			

		}

		//printf("I'm about to enter the sec loop \n");
		//looking for pyx, wyx and posw
		//DEBUGG
                for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
                //for(int x = 123+spat_patch_step; x <123+spat_patch_step+1 ; x++)

		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			//DEBUGG
                        for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
                        //for(int y = 118+spat_patch_step; y < 118+spat_patch_step+1; y++)
			{

				//linear index of the current pixel at the center of the search window
				int l_2d=(y-spat_patch_step)+lx;

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//auxiliary vector
				std::vector<float> tmp_wyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));

				//(re)initializing variables
				int w=0;
			
				//DEBUGG
				//printf("I'm here now \n");

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					//DEBUGG
					//printf("yo 1\n");

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//the linear index of the current pixel in the search window
						int l_ser=(y_ser-spat_patch_step)+lx_ser+lh;
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						 //DEBUGG
			                         /*printf("here I am \n");
						 std::cout<<"x= "<<x<<", y= "<<y<<", x_ser= "<<x_ser<<", y_ser= "<<y_ser<<", lx_ser= "<<lx_ser<<", l_ser= "<<l_ser<<", l_ser_2d= "<<l_ser_2d<<std::endl;
						 std::cout<<"l_ser= "<<l_ser<<", l_2d= "<<l_2d<<", l_ser_2d= "<<l_ser_2d<<std::endl;*/

						/*for(int i=0; i<pxy[l_ser].size();i++)
						{
							std::cout<<"pxy= "<<pxy[l_ser][i]<<std::endl;
						}*/

						//we verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
						int pos=std::find(pxy[l_ser].begin(), pxy[l_ser].end(),l_2d)-pxy[l_ser].begin();

						//DEBUGG
						//printf("things found \n");

						//in the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
						if( pos < (int) pxy[l_ser].size() )
						{
							tmp_wyx[w]=wxy[l_ser][pos];
							tmp_pyx[w]=l_ser_2d;
							tmp_posw[w]=pos;
							w++;

							//DEBUGG
							//printf("coucou 1\n");
						}
					 //DEBUGG
                                         //printf("coucou 2\n");

					}
					//DEBUGG
                                        //printf("coucou 3\n");


				}

				//DEBUGG
				//printf("about to store things\n");
				tmp_wyx.resize(w);
				tmp_pyx.resize(w);
				tmp_posw.resize(w);

				wyx.push_back(tmp_wyx);
				pyx.push_back(tmp_pyx);
				posw.push_back(tmp_posw);

				//DEBUGG
				//printf("stored \n");

			}

		}



	}

	//printf("If I'm here it means that I'm going out soone \n");

}





//This routines computes 3D distances on the multispectral image weighted with coefficient of the Smatrix
//the distances are unclassified
void nlweights_3Dall_ms_Smatrix_unclassified(float *ms_image_new_ptr, std::vector<float> &vect_xy, float hSim, float hClose,
		int reswind, int spat_patch_step, int ms_channels,int width_ms_new,
		int height_ms_new,float **Smatrix,int x, int y, int num_hs_channel)
{

	// width_ms_new and height_ms_new are the spatial dimensions of the multispectral image
	// after spatial extension


	// Adapt filter parameters to size of comparison window
	int patch_size_3d =(2 * spat_patch_step + 1)*(2*spat_patch_step+1)*(ms_channels);//the size of the 3d patch

	float hSim2_3d = hSim * hSim * patch_size_3d;

	float hClose2 = hClose * hClose;

	int dim=(2*spat_patch_step + 1)*(2*spat_patch_step+1);

	std::vector<float>  diff_3d(patch_size_3d);
	std::vector<float>  diff_2_3d(patch_size_3d);

	std::vector<float>  tmp_pix_central_3d(patch_size_3d);

	//Restraining the search in the spatial (x) domain in case we are near the borders
	int bound_start_x=MAX(x-reswind,spat_patch_step);
	int bound_end_x=MIN(x+reswind,height_ms_new-spat_patch_step-1);

	//Restraining the search in the spatial (y) domain in case we are near the borders
	int bound_start_y=MAX(y-reswind,spat_patch_step);
	int bound_end_y=MIN(y+reswind,width_ms_new-spat_patch_step-1);


	//slicing (isolating) the 3d patch out of the image
	float *coef_tab=new float[ms_channels];
	float sum_coef=0.0;
	for(int m=0;m<ms_channels;m++)
	{
		float sval=Smatrix[m][num_hs_channel];
		coef_tab[m]=sval;
		sum_coef=sum_coef+sval;
	}

	int ind_cent=0;
	for(int m=0;m<ms_channels;m++)
	{
		float coef=coef_tab[m];
		array_slicing_coef(x-spat_patch_step,y-spat_patch_step,m,x+spat_patch_step
				,y+spat_patch_step,m,height_ms_new,width_ms_new,ms_image_new_ptr,ind_cent,tmp_pix_central_3d,coef);
		ind_cent=ind_cent+dim;
	}

	//(re)initializing variables
	int w=0; //browses the search window for each central pixel
	int wcentral=0;
	float dMin_3d = fLarge;
	float dSim_3d=0;
	float dClose=0;

	//entering the search window
	for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
	{

		for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
		{

			std::vector<float>  tmp_pix_ser_3d(patch_size_3d);

			//slicing (isolating) the 3d patch out of the image
			int ind_ser=0;
			for(int m=0;m<ms_channels;m++)
			{
				float coef=coef_tab[m];
				array_slicing_coef(x_ser-spat_patch_step,y_ser-spat_patch_step,m,x_ser+spat_patch_step
						,y_ser+spat_patch_step,m,height_ms_new,width_ms_new,ms_image_new_ptr,ind_ser,tmp_pix_ser_3d,coef);
				ind_ser=ind_ser+dim;
			}

			//BILATERAL DISTANCE WITH THE 3D PATCHS
			//computing the patch-based similarity distance
			substract_two_vectors(tmp_pix_central_3d,tmp_pix_ser_3d,patch_size_3d,diff_3d);

			power_two_vector(diff_3d,patch_size_3d,diff_2_3d);

			dSim_3d=sum_values_vector(diff_2_3d,patch_size_3d);

			float hSim2=hSim2_3d*sum_coef;
			dSim_3d /= hSim2;

			//computing the spatial distance
			dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
			dClose /= hClose2;


			//storing the bilateral distances computed with the 3d patches and 2d patches
			float d_bilateral_3d=dSim_3d+dClose;

			vect_xy[w]=d_bilateral_3d;

			/// Save w-position of central pixel in neighborhood and the minimum bilateral distance
			if((x==x_ser) && (y==y_ser))
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

	//restraining the size of WandPxy[l] to w
	vect_xy.resize(w);

	//Assign minimum distance to central pixel
	vect_xy[wcentral]=dMin_3d;
}




// Computes the sqrt(weights) on the full upsampled HS image, so 3D patches are considered.
// The patch distance across channels is coupled with the L2 norm.
// The same weight distribution is obtained for all HS channels.
// Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
void nlweights_hs_3Dall(libUSTG::cflimage hyperDataUpsampled, std::vector< std::vector<float> > &wxy,
		std::vector< std::vector<int> > &posxy, std::vector< std::vector<float> > &wyx,
		std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw, float hSim, float hClose,
		int numSimPixels, int reswind, int spat_patch_step, int flagNormalize)
{

	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	//extension of the image
	//in this case the image should be extended only spatially
	libUSTG::cflimage im_3d_new=extend_image_dimensions_spat(hyperDataUpsampled,spat_patch_step); //spatial extension
	float *im_3d_new_ptr=im_3d_new.d_v;

	//recovering the dimensions of the images
	int width=hyperDataUpsampled.w();
	int hs_channels=hyperDataUpsampled.c();

	int width_new=im_3d_new.w();
	int height_new=im_3d_new.h();

	// Adapt filter parameters to size of comparison window
	int patch_size_3d =(2 * spat_patch_step + 1)*(2*spat_patch_step+1)*(hs_channels);//the size of the 3d patch

	float hSim2_3d = hSim * hSim * patch_size_3d;

	float hClose2 = hClose * hClose;

	std::vector<float>  diff_3d(patch_size_3d);
	std::vector<float>  diff_2_3d(patch_size_3d);

	std::vector<float>  tmp_pix_central_3d(patch_size_3d);
	std::vector<float>  tmp_pix_ser_3d(patch_size_3d);


	//browsing the rows
	//the height here is the NUMBER OF ROWS

	for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
	{
		int lx=(x-spat_patch_step)*width;

		//Restraining the search in the spatial (x) domain in case we are near the borders
		int bound_start_x=MAX(x-reswind,spat_patch_step);
		int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

		//browsing the columns
		//the width here is THE NUMBER OF COLUMNS
		for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
		{

			int l_2d=lx+(y-spat_patch_step);

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_y=MAX(y-reswind,spat_patch_step);
			int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

			//slicing (isolating) the 3d patch out of the image
			array_slicing(x-spat_patch_step,y-spat_patch_step,0,x+spat_patch_step
					,y+spat_patch_step,hs_channels-1,height_new,width_new,im_3d_new_ptr,tmp_pix_central_3d);

			//(re)initializing variables
			int w=0; //browses the search window for each central pixel
			int wcentral=0;
			float dMin_3d = fLarge;
			float dSim_3d=0;
			float dClose=0;

			//auxiliary vector
			std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));
			std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));

			//entering the search window
			for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
			{
				int lx_ser = (x_ser-spat_patch_step) * width;

				for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
				{

					//linear index of the current pixel in the search window
					int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

					//slicing (isolating) the 3d patch out of the image
					array_slicing(x_ser-spat_patch_step,y_ser-spat_patch_step,0,x_ser+spat_patch_step,
							y_ser+spat_patch_step,hs_channels-1,height_new,width_new,im_3d_new_ptr,tmp_pix_ser_3d);

					//BILATERAL DISTANCE WITH THE 3D PATCHS
					//computing the patch-based similarity distance
					substract_two_vectors(tmp_pix_central_3d,tmp_pix_ser_3d,patch_size_3d,diff_3d);

					power_two_vector(diff_3d,patch_size_3d,diff_2_3d);

					dSim_3d=sum_values_vector(diff_2_3d,patch_size_3d);

					dSim_3d /= hSim2_3d;

					//computing the spatial distance
					dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
					dClose /= hClose2;

					//storing the bilateral distances computed with the 3d patches and 2d patches
					float d_bilateral_3d=dSim_3d+dClose;
					tmp_wxy[w]=d_bilateral_3d;

					// Save w-position of central pixel in neighborhood and the minimum bilateral distance
					if((x==x_ser) && (y==y_ser))
					{
						wcentral = w;
					}
					else
					{
						if(d_bilateral_3d < dMin_3d)
							dMin_3d = d_bilateral_3d;

					}

					// Save neighbouring pixel position
					tmp_pxy[w] = l_ser_2d;

					// Update index
					w++;


				}
			}

			//restraining the size of tmp_wxy and tmp_pxy to w
			tmp_wxy.resize(w);
			tmp_pxy.resize(w);

			//Assign minimum distance to central pixel
			tmp_wxy[wcentral]=dMin_3d;

			//Adapt number of neighboring pixels to window size
			int numSimPixels0 = MIN(numSimPixels, w);

			// Order bilateral distances with the 3d patch only
			sort_vects(tmp_wxy, tmp_pxy);

			//Keep only the first numSimPixels0 from wxy and posxy
			tmp_wxy.resize(numSimPixels0);
			tmp_pxy.resize(numSimPixels0);

			// applying the exponential
			float sum_weights=0.0f;
			for(int r = 0; r < numSimPixels0; r++)
			{
				float weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
				tmp_wxy[r]=weight;
				sum_weights += weight;
			}

			//normalization and applying the square root (sqrt)
			if(flagNormalize)
			{
				if(sum_weights > fTiny)
				{
					for(int r = 0; r < numSimPixels0; r++)
					{
						float weight = (float) tmp_wxy[r] / sum_weights;
						tmp_wxy[r] = sqrt(weight);
					}

				} else
				{
					tmp_wxy.clear();
					tmp_pxy.clear();

					tmp_wxy.push_back(1.0f);
					tmp_pxy.push_back(l_2d);
				}

			} else
			{
				for(int r = 0; r < numSimPixels0; r++)
				{
					float weight = (float)tmp_wxy[r];
					tmp_wxy[r] = sqrt(weight);
				}
			}

			wxy.push_back(tmp_wxy);
			posxy.push_back(tmp_pxy);
		}
	}


	//looking for posxy, wyx and posw
	for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
	{
		int lx = (x-spat_patch_step) * width;

		//Restraining the search in the spatial (x) domain in case we are near the borders
		int bound_start_x=MAX(x-reswind,spat_patch_step);
		int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

		//browsing the columns
		//the width here is THE NUMBER OF COLUMNS
		for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
		{

			//linear index of the current pixel at the center of the search window
			int l_2d=(y-spat_patch_step)+lx;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_y=MAX(y-reswind,spat_patch_step);
			int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

			//auxiliary vector
			std::vector<float> tmp_wyx((2*reswind+1)*(2*reswind+1));
			std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
			std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));

			//(re)initializing variables
			int w=0;
			//entering the search window
			for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
			{
				int lx_ser = (x_ser-spat_patch_step) * width;

				for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
				{

					//the linear index of the current pixel in the search window
					int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

					//we verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
					int pos=std::find(posxy[l_ser_2d].begin(), posxy[l_ser_2d].end(),l_2d)-posxy[l_ser_2d].begin();

					//in the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
					//WeightsAndPositionsYX s;
					if( pos < (int) posxy[l_ser_2d].size() )
					{
						tmp_wyx[w]=wxy[l_ser_2d][pos];
						tmp_pyx[w]=l_ser_2d;
						tmp_posw[w]=pos;
						w++;

					}

				}
			}
			tmp_wyx.resize(w);
			tmp_pyx.resize(w);
			tmp_posw.resize(w);

			wyx.push_back(tmp_wyx);
			posyx.push_back(tmp_pyx);
			posw.push_back(tmp_posw);
		}

	}

}



// Computes the sqrt(weights) using the nearest bands w.r.t. the central band, so 3D patches are considered.
// The patch distance across channels is coupled with the L2 norm.
// A different weight distribution is obtained for each HS channel.
// Weight structure: wxy[h*dim+i][n] for each channel h, each pixel i, and each neighbour n.
void nl_weights_hs_3Dnearest(libUSTG::cflimage hyperDataUpsampled, std::vector< std::vector<float> > &wxy,
		std::vector< std::vector<int> > &posxy, std::vector< std::vector<float> > &wyx,
		std::vector< std::vector<int> > &posyx, std::vector< std::vector<int> > &posw, float hSim, float hClose,
		int numSimPixels, int reswind, int spat_patch_step, int spec_patch_step, int flagNormalize)
{

	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	// Extension of the image
	libUSTG::cflimage im_3d_new_spat = extend_image_dimensions_spat(hyperDataUpsampled, spat_patch_step); //spatial extension
	libUSTG::cflimage im_3d_new = extend_image_dimensions_spec(im_3d_new_spat, spec_patch_step); //spectral extension

	// Pointer to the values of the extended image
	float *im_3d_new_ptr = im_3d_new.d_v;

	// Recovering the dimensions of the images
	int width = hyperDataUpsampled.w();
	int height = hyperDataUpsampled.h();
	int width_new = im_3d_new.w();
	int height_new = im_3d_new.h();
	int hs_channels_new = im_3d_new.c();

	// Adapt filter parameters to size of comparison window
	int patch_size_3d = (2 * spat_patch_step + 1) * (2 * spat_patch_step + 1) * (2 * spec_patch_step + 1); // Size of the 3D patch
	float hSim2_3d = hSim * hSim * patch_size_3d;
	float hClose2 = hClose * hClose;


	// Auxiliary vectors for the computations
	std::vector<float> diff_3d(patch_size_3d);
	std::vector<float> diff_2_3d(patch_size_3d);
	std::vector<float> tmp_pix_central_3d(patch_size_3d);
	std::vector<float> tmp_pix_ser_3d(patch_size_3d);


	for(int h = spec_patch_step; h < hs_channels_new-spec_patch_step; h++)
	{
		int lh = (h-spec_patch_step)*width*height;

		//browsing the rows
		//the height here is the NUMBER OF ROWS
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		{

			int lx=(x-spat_patch_step)*width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			{

				//int l=lh+lx+(y-spat_patch_step);
				int l_2d=lx+(y-spat_patch_step);

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//slicing (isolating) the 3d patch (centered on the central pixel (x,y,h)) out of the image
				array_slicing(x-spat_patch_step,y-spat_patch_step,h-spec_patch_step,x+spat_patch_step
						,y+spat_patch_step,h+spec_patch_step,height_new,width_new,im_3d_new_ptr,tmp_pix_central_3d);


				//(re)initializing variables
				int w=0; //browses the search window for each central pixel
				int wcentral=0;
				float dMin_3d = fLarge;
				float dSim_3d=0;
				float dClose=0;

				//auxiliary vector
				std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//linear index of the current pixel in the search window
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						//slicing (isolating) the 3d patch (centered on the pixel (x_ser,y_ser,h)) out of the image
						array_slicing(x_ser-spat_patch_step,y_ser-spat_patch_step,h-spec_patch_step,x_ser+spat_patch_step,
								y_ser+spat_patch_step,h+spec_patch_step,height_new,width_new,im_3d_new_ptr,tmp_pix_ser_3d);


						//BILATERAL DISTANCE WITH THE 3D PATCHS
						//computing the patch-based similarity distance
						substract_two_vectors(tmp_pix_central_3d,tmp_pix_ser_3d,patch_size_3d,diff_3d);

						power_two_vector(diff_3d,patch_size_3d,diff_2_3d);

						dSim_3d=sum_values_vector(diff_2_3d,patch_size_3d);

						dSim_3d /= hSim2_3d;

						//computing the spatial distance
						dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
						dClose /= hClose2;


						//storing the bilateral distances computed with the 3d patches and 2d patches
						float d_bilateral_3d=dSim_3d+dClose;
						tmp_wxy[w]=d_bilateral_3d;

						// Save w-position of the central pixel in neighborhood and the minimum bilateral distance
						if((x==x_ser) && (y==y_ser))
						{
							wcentral = w;
						}
						else
						{
							if(d_bilateral_3d < dMin_3d)
								dMin_3d = d_bilateral_3d;
						}

						//Save neighbouring pixel position
						tmp_pxy[w]=l_ser_2d;

						// Update index
						w++;


					}
				}

				//restraining the size of tmp_wxy and tmp_pxy to w
				tmp_wxy.resize(w);
				tmp_pxy.resize(w);

				//Assign minimum distance to central pixel
				tmp_wxy[wcentral]=dMin_3d;

				//Adapt number of neighboring pixels to window size
				int numSimPixels0 = MIN(numSimPixels, w);

				//sorting the values and there positions wrt the ascending order of the value
				sort_vects(tmp_wxy, tmp_pxy);

				//Keep only the first numSimPixels0 from wxy and posxy
				tmp_wxy.resize(numSimPixels0);
				tmp_pxy.resize(numSimPixels0);

				// applying the exponential
				float sum_weights=0.0f;
				for(int r = 0; r < numSimPixels0; r++)
				{
					float weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
					tmp_wxy[r]=weight;
					sum_weights += weight;
				}

				//normalization and applying the square root (sqrt)
				if(flagNormalize)
				{
					if(sum_weights > fTiny)
					{
						for(int r = 0; r < numSimPixels0; r++)
						{
							float weight = (float) tmp_wxy[r] / sum_weights;
							tmp_wxy[r] = sqrt(weight);
						}

					} else
					{
						tmp_wxy.clear();
						tmp_pxy.clear();

						tmp_wxy.push_back(1.0f);
						tmp_pxy.push_back(l_2d);
					}

				}
				else
				{
					for(int r = 0; r < numSimPixels0; r++)
					{
						float weight = (float)tmp_wxy[r];
						tmp_wxy[r] = sqrt(weight);
					}
				}

				wxy.push_back(tmp_wxy);
				posxy.push_back(tmp_pxy);

			}
		}


		//looking for posyx, wyx and posw
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			{

				//linear index of the current pixel at the center of the search window
				//int l=(y-spat_patch_step)+lx+lh;
				int l_2d=(y-spat_patch_step)+lx;

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//auxiliary vector
				std::vector<float> tmp_wyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));

				//(re)initializing variables
				int w=0;

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//the linear index of the current pixel in the search window
						int l_ser=(y_ser-spat_patch_step)+lx_ser+lh;
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						//we verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
						int pos = std::find(posxy[l_ser].begin(), posxy[l_ser].end(),l_2d) - posxy[l_ser].begin();

						//in the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
						if(pos < (int) posxy[l_ser].size())
						{
							tmp_wyx[w]=wxy[l_ser][pos];
							tmp_pyx[w]=l_ser_2d;
							tmp_posw[w]=pos;
							w++;

						}

					}
				}
				tmp_wyx.resize(w);
				tmp_pyx.resize(w);
				tmp_posw.resize(w);

				wyx.push_back(tmp_wyx);
				posyx.push_back(tmp_pyx);
				posw.push_back(tmp_posw);

			}

		}


	}

}




// Computes the sqrt(weights) by first choosing most similar neighbours on the full upsampled HS image (3D patches) and computing then
// the weights using 2D patches centered at the pre-selected pixels on the corresponding band.
// A different weight distribution is obtained for each HS channel.
// Weight structure: wxy[h*dim+i][n] for each channel h, each pixel i, and each neighbour n.
void nlweights_hs_3Dclassify_2Dband(libUSTG::cflimage hyperDataUpsampled, std::vector< std::vector<float> > &wxy,
		std::vector< std::vector<int> > &pxy, std::vector< std::vector<float> > &wyx,
		std::vector< std::vector<int> > &pyx, std::vector< std::vector<int> > &posw, float hSim,
		float hClose, int numSimPixels, int reswind, int spat_patch_step, int flagNormalize)
{


	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	//extension of the image
	//in this case the image should be extended only spatially
	libUSTG::cflimage im_3d_new=extend_image_dimensions_spat(hyperDataUpsampled,spat_patch_step); //spatial extension
	float *im_3d_new_ptr=im_3d_new.d_v;

	//recovering the dimensions of the images
	int width=hyperDataUpsampled.w();
	int height=hyperDataUpsampled.h();
	int hs_channels=hyperDataUpsampled.c();

	int width_new=im_3d_new.w();
	int height_new=im_3d_new.h();


	// Adapt filter parameters to size of comparison window
	int patch_size_2d =(2 * spat_patch_step + 1)*(2*spat_patch_step+1);//the size of the 2d patch

	float hSim2_2d = hSim * hSim * patch_size_2d;

	float hClose2 = hClose * hClose;

	std::vector<float>  diff_2d(patch_size_2d);
	std::vector<float>  diff_2_2d(patch_size_2d);

	std::vector<float>  tmp_pix_central_2d(patch_size_2d);
	std::vector<float>  tmp_pix_ser_2d(patch_size_2d);

	//declaring the matrix where 3d unclassified weights computed with all the spectral bands are going to be stored
	std::vector< std::vector<float> > wxy_3D;

	//Computing the unclassified 3d weights
	nlweights_3Dall_unclassified(im_3d_new_ptr, wxy_3D, hSim, hClose, reswind, spat_patch_step, hs_channels, width_new,
			height_new);

	for(int h=0;h<hs_channels;h++)
	{
		int lh = h*width*height;
		//browsing the rows
		//the height here is the NUMBER OF ROWS
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			{
				//linear index of the current pixel at the center of the search window
				int l_comp=(y-spat_patch_step)+lx;//index for the comparison matrix

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//slicing (isolating) the 2d patch (centered on the central pixel (x,y,0)) out of the image
				array_slicing(x-spat_patch_step,y-spat_patch_step,h,x+spat_patch_step
						,y+spat_patch_step,h,height_new,width_new,im_3d_new_ptr,tmp_pix_central_2d);

				//(re)initializing variables
				int w=0; //browses the search window for each central pixel
				int wcentral=0;
				float dMin_2d = fLarge;
				float dSim_2d=0;
				float dClose=0;

				//auxiliary vector
				std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pxy((2*reswind+1)*(2*reswind+1));

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//linear index of the current pixel in the search window
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;


						//slicing (isolating) the 2d patch out of the image
						array_slicing(x_ser-spat_patch_step,y_ser-spat_patch_step,h,x_ser+spat_patch_step,
								y_ser+spat_patch_step,h,height_new,width_new,im_3d_new_ptr,tmp_pix_ser_2d);

						//BILATERAL DISTANCE WITH THE 2D PATCHS
						//computing the patch-based similarity distance
						substract_two_vectors(tmp_pix_central_2d,tmp_pix_ser_2d,patch_size_2d,diff_2d);

						power_two_vector(diff_2d,patch_size_2d,diff_2_2d);

						dSim_2d=sum_values_vector(diff_2_2d,patch_size_2d);

						dSim_2d /= hSim2_2d;

						//computing the spatial distance
						dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
						dClose /= hClose2;


						//storing the bilateral distances computed with the 3d patches and 2d patches
						float d_bilateral_2d=dSim_2d+dClose;
						tmp_wxy[w]=d_bilateral_2d;

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

						//Save neighbouring pixel position
						tmp_pxy[w]=l_ser_2d;

						// Update index
						w++;


					}
				}

				//restraining the size of tmp_wxy to w
				tmp_wxy.resize(w);
				tmp_pxy.resize(w);


				//Assign minimum distance to central pixel
				tmp_wxy[wcentral]=dMin_2d;

				//Adapt number of neighboring pixels to window size
				int numSimPixels0 = MIN(numSimPixels, w);

				std::vector<float> vec=wxy_3D[l_comp];
				sort_wrt_comp_vect(tmp_wxy, tmp_pxy, vec);

				//Keep only the first numSimPixels0 from wxy and posxy
				tmp_wxy.resize(numSimPixels0);
				tmp_pxy.resize(numSimPixels0);

				// applying the exponential
				float sum_weights=0.0f;
				for(int r = 0; r < numSimPixels0; r++)
				{
					float weight = libUSTG::wxSLUT(tmp_wxy[r], lut);
					tmp_wxy[r]=weight;
					sum_weights += weight;
				}

				//normalization and applying the square root (sqrt)
				if(flagNormalize)
				{
					if(sum_weights > fTiny)
					{
						for(int r = 0; r < numSimPixels0; r++)
						{
							float weight = (float) tmp_wxy[r] / sum_weights;
							tmp_wxy[r] = sqrt(weight);
						}

					} else
					{
						tmp_wxy.clear();
						tmp_pxy.clear();

						tmp_wxy.push_back(1.0f);
						tmp_pxy.push_back(l_comp);
					}

				} else
				{
					for(int r = 0; r < numSimPixels0; r++)
					{
						float weight = (float) tmp_wxy[r];
						tmp_wxy[r] = sqrt(weight);
					}
				}

				wxy.push_back(tmp_wxy);
				pxy.push_back(tmp_pxy);
			}
		}

		//looking for pyx, wyx and posw
		for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
		{
			int lx = (x-spat_patch_step) * width;

			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_x=MAX(x-reswind,spat_patch_step);
			int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

			//browsing the columns
			//the width here is THE NUMBER OF COLUMNS
			for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
			{

				//linear index of the current pixel at the center of the search window
				//int l=(y-spat_patch_step)+lx+lh;
				int l_2d=(y-spat_patch_step)+lx;

				//Restraining the search in the spatial (x) domain in case we are near the borders
				int bound_start_y=MAX(y-reswind,spat_patch_step);
				int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

				//auxiliary vector
				std::vector<float> tmp_wyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_pyx((2*reswind+1)*(2*reswind+1));
				std::vector<int> tmp_posw ((2*reswind+1)*(2*reswind+1));

				//(re)initializing variables
				int w=0;

				//entering the search window
				for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
				{
					int lx_ser = (x_ser-spat_patch_step) * width;

					for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
					{

						//the linear index of the current pixel in the search window
						int l_ser=(y_ser-spat_patch_step)+lx_ser+lh;
						int l_ser_2d=(y_ser-spat_patch_step)+lx_ser;

						//we verify if the pixels in the search window have (x,y) the current pixel in the search window as a neighbour
						int pos=std::find(pxy[l_ser].begin(), pxy[l_ser].end(),l_2d)-pxy[l_ser].begin();

						//in the case where (x,y) might be a possible neighbor to (x_ser,y_ser)
						if( pos < (int) pxy[l_ser].size() )
						{
							tmp_wyx[w]=wxy[l_ser][pos];
							tmp_pyx[w]=l_ser_2d;
							tmp_posw[w]=pos;
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





// Computes weights on the full upsampled HS image, so 3D patches are considered, but without classifiying them.
// The same weight is obtained for all HS channels.
// Weight structure: wxy[i][n] for each pixel i, and each neighbour n.
void nlweights_3Dall_unclassified(float *im_3d_new_ptr, std::vector< std::vector<float> > &wxy, float hSim, float hClose,int reswind,
		int spat_patch_step, int hs_channels, int width_new, int height_new)
{

	// Tabulate the function Exp(-x) for x>0.
	int luttaille = (int) (LUTMAX * LUTPRECISION);
	float *lut = new float[luttaille];
	libUSTG::wxFillExpLut(lut, luttaille);

	// Adapt filter parameters to size of comparison window
	int patch_size_3d =(2 * spat_patch_step + 1)*(2*spat_patch_step+1)*(hs_channels);//the size of the 3d patch

	float hSim2_3d = hSim * hSim * patch_size_3d;

	float hClose2 = hClose * hClose;

	std::vector<float>  diff_3d(patch_size_3d);
	std::vector<float>  diff_2_3d(patch_size_3d);

	std::vector<float>  tmp_pix_central_3d(patch_size_3d);
	std::vector<float>  tmp_pix_ser_3d(patch_size_3d);


	/*The 3d weights are computed with all the spectral bands. Since these 3d weights are going to be the same
	 *in all bands for each spatial coordinates, these computations are done once
	 */

	//browsing the rows
	//the height here is the NUMBER OF ROWS
	for(int x = spat_patch_step; x < height_new-spat_patch_step; x++)
	{


		//Restraining the search in the spatial (x) domain in case we are near the borders
		int bound_start_x=MAX(x-reswind,spat_patch_step);
		int bound_end_x=MIN(x+reswind,height_new-spat_patch_step-1);

		//browsing the columns
		//the width here is THE NUMBER OF COLUMNS
		for(int y = spat_patch_step; y < width_new-spat_patch_step; y++)
		{
			//Restraining the search in the spatial (x) domain in case we are near the borders
			int bound_start_y=MAX(y-reswind,spat_patch_step);
			int bound_end_y=MIN(y+reswind,width_new-spat_patch_step-1);

			//slicing (isolating) the 3d patch out of the image
			array_slicing(x-spat_patch_step,y-spat_patch_step,0,x+spat_patch_step
					,y+spat_patch_step,hs_channels-1,height_new,width_new,im_3d_new_ptr,tmp_pix_central_3d);

			//(re)initializing variables
			int w=0; //browses the search window for each central pixel
			int wcentral=0;
			float dMin_3d = fLarge;
			float dSim_3d=0;
			float dClose=0;

			//auxiliary vector
			std::vector<float> tmp_wxy((2*reswind+1)*(2*reswind+1));

			//entering the search window
			for(int x_ser=bound_start_x; x_ser<= bound_end_x; x_ser++) //x_ser (x in search window)
			{

				for(int y_ser=bound_start_y; y_ser <= bound_end_y; y_ser++) //y_ser (y in search window)
				{

					//slicing (isolating) the 3d patch out of the image
					array_slicing(x_ser-spat_patch_step,y_ser-spat_patch_step,0,x_ser+spat_patch_step,
							y_ser+spat_patch_step,hs_channels-1,height_new,width_new,im_3d_new_ptr,tmp_pix_ser_3d);

					//BILATERAL DISTANCE WITH THE 3D PATCHS
					//computing the patch-based similarity distance
					substract_two_vectors(tmp_pix_central_3d,tmp_pix_ser_3d,patch_size_3d,diff_3d);

					power_two_vector(diff_3d,patch_size_3d,diff_2_3d);

					dSim_3d=sum_values_vector(diff_2_3d,patch_size_3d);

					dSim_3d /= hSim2_3d;

					//computing the spatial distance
					dClose=(float) ((x - x_ser) * (x - x_ser) + (y - y_ser) * (y - y_ser));
					dClose /= hClose2;


					//storing the bilateral distances computed with the 3d patches and 2d patches
					float d_bilateral_3d=dSim_3d+dClose;

					tmp_wxy[w]=d_bilateral_3d;

					// Save w-position of central pixel in neighborhood and the minimum bilateral distance
					if((x==x_ser) && (y==y_ser))
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

			//restraining the size of tmp_wxy to w
			tmp_wxy.resize(w);

			//Assign minimum distance to central pixel
			tmp_wxy[wcentral]=dMin_3d;

			wxy.push_back(tmp_wxy);

		}
	}
}





// Sorts two vectors wrt the values vector (ascending order)
void sort_vects(std::vector<float> &values_vect, std::vector<int> &positions_vect)
{
	int size=values_vect.size();
	std::vector<ValAndPos> vector(size);

	// Copy everything in a structure
	for(int i=0;i<size;i++)
	{
		vector[i].value=values_vect[i];
		vector[i].position=positions_vect[i];
	}

	// Sorting (in the ascending direction)
	std::sort(vector.begin(), vector.end(), [](const ValAndPos &a, const ValAndPos &b) { return a.value < b.value; } );

	//printf("nearly out \n");
	//recopy everything into the original vectors
	for(int i=0;i<size;i++)
	{
		values_vect[i]=vector[i].value;
		positions_vect[i]=vector[i].position;
	}

}



// Sort with respect to a comparison vector (ascending order)
void sort_wrt_comp_vect(std::vector<float> &values_vect, std::vector<int> &positions_vect, const std::vector<float> &comparison_vect)
{
	int size=values_vect.size();
	std::vector<ValPosCompVect> vector(size);

	//copy everything in a structure
	for(int i=0;i<size;i++)
	{
		vector[i].value=values_vect[i];
		vector[i].comp_value=comparison_vect[i];
		vector[i].position=positions_vect[i];
	}

	//sorting(in the ascending direction)
	std::sort(vector.begin(), vector.end(),
			[](const ValPosCompVect &a, const ValPosCompVect &b) { return a.comp_value < b.comp_value; } );

	//recopy everything into the original vectors
	for(int i=0;i<size;i++)
	{
		values_vect[i]=vector[i].value;
		//comparison_vect[i]=vector[i].comp_value;
		positions_vect[i]=vector[i].position;
	}

}

void array_slicing_coef(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height, int width,
		float *im_to_slice, int ind, std::vector<float>  &sliced_vector,float coef)
{
	int cpt=ind;
	for(int h=h_start;h<=h_end;h++)
	{
		for(int i=i_start;i<=i_end;i++)
		{
			for(int j=j_start;j<=j_end;j++)
			{
				sliced_vector[cpt]=coef*im_to_slice[j+(i*width)+(h*height*width)];
				cpt++;
			}
		}
	}
}



void array_slicing(int i_start, int j_start, int h_start, int i_end, int j_end, int h_end, int height, int width,
		float *im_to_slice, std::vector<float>  &sliced_vector)
{

	int cpt=0;

	for(int h=h_start;h<=h_end;h++)
	{
		for(int i=i_start;i<=i_end;i++)
		{
			for(int j=j_start;j<=j_end;j++)
			{
				sliced_vector[cpt]=im_to_slice[j+(i*width)+(h*height*width)];
				cpt++;
			}
		}
	}
}



float sum_values_vector(std::vector<float> &vec, int size_vec)
{
	float sum=0.0;

	for(int j=0;j<size_vec;j++)
	{
		sum=sum+vec[j];
	}

	return sum;
}



void substract_two_vectors(std::vector<float> &vec1, std::vector<float> &vec2, int size_vec, std::vector<float> &res_vector)
{

	for(int i=0;i<size_vec;i++)
	{
		res_vector[i]=vec1[i]-vec2[i];
	}

}



void power_two_vector(std::vector<float> &vector, int size_vec, std::vector<float> &res_vector)
{
	for(int k=0;k<size_vec;k++)
	{
		//res_vector[k]=pow(vector[k],2);
		res_vector[k]=vector[k]*vector[k];
	}
}



// Prints the values of a 1D pointer array.
void print_values_array_1D_pt(float *array, int size)
{

	for(int i=0;i<size;i++)
	{

		printf(" %.4f",array[i]);
	}

	printf("\n");
}


// Prints the values of an array from the cflimage class.
void print_values_array(libUSTG::cflimage vec)
{
	int nr=vec.h(); //number of rows (height)
	int nc=vec.w(); //number of columns (width)
	int nb=vec.c(); //number of channels

	//
	for(int h=0;h<nb;h++)
	{
		for(int i=0;i<nr;i++)
		{
			for(int j=0;j<nc;j++)
			{
				printf(" %.4f",vec[j+(i*nc)+(h*nc*nr)]);
			}
			printf("\n");
		}
		printf("\n");
	}
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



// Extends the SPECTRAL dimensions of an image by mirror effect which helps in case the treatment of the image takes the borders into
// account. The image to extend should be a 3d image.
libUSTG::cflimage  extend_image_dimensions_spec(libUSTG::cflimage im_to_extend, int spec_patch_step)
{

	/*extracting the dimensions of the image to extend
	 *b_old is the number of the bands of the image to extend
	 *The new b is the one of the extended image
	 */

	int nr=im_to_extend.h(); //number of rows of the image to extend (it's not going to change)
	int nc=im_to_extend.w(); //number of columns of the image to extend (it's not going to change)
	int b_old=im_to_extend.c(); //number of spectral bands of the image to extend

	/*Test in case the patch step is bigger than the number of bands of the image to extend
	 *If it's the case (which I don't find interesting), it's considered as an error and the
	 *image is not extended.
	 */

	if((spec_patch_step>b_old))
	{
		printf("The patch step is bigger than the number of bands if the image to extend, not an interesting case \n");
		return im_to_extend;
	}


	/*Test in case the patch step is null*/
	if(spec_patch_step==0)
	{
		return im_to_extend;
	}

	/* step1: creating a matrix with the new spectral dimension */
	int b_new=b_old+(2*spec_patch_step);

	libUSTG::cflimage extended_image;
	extended_image.create(nc,nr,b_new);

	int k,j,i,m=0;
	/*step2: copy of the content of the old matrix into the new matrix*/
	/*the browsing is done in rows*/

	for(k=spec_patch_step;k<=b_new-spec_patch_step-1;k++)
	{
		for(j=0;j<nc;j++)
		{
			for(i=0;i<nr;i++)
			{
				extended_image[(k*nc*nr)+(j*nr)+i]=im_to_extend[m];
				m++;
			}
		}

	}

	/*filling the front bands*/
	int k_prime=(2*spec_patch_step)-1;
	for(k=0;k<=spec_patch_step-1;k++)
	{
		for(j=0;j<nc;j++)
		{
			for(i=0;i<nr;i++)
			{
				extended_image[(k*nc*nr)+(j*nr)+i]=extended_image[(k_prime*nc*nr)+(j*nr)+i];
			}
		}

		k_prime--;
	}

	/*filling the back rows*/
	int k_prime_second=b_new-(2*spec_patch_step);
	for(k=b_new-1;k>=b_new-spec_patch_step;k--)
	{
		for(j=0;j<nc;j++)
		{
			for(i=0;i<nr;i++)
			{
				extended_image[(k*nc*nr)+(j*nr)+i]=extended_image[(k_prime_second*nc*nr)+(j*nr)+i];
			}
		}

		k_prime_second++;

	}

	return extended_image;
}

} // libUSTGHYPER
