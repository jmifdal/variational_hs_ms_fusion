hdrdir=   -I/usr/local/include/
libdir=    -L/usr/local/lib/     

COPT	= -O3 -funroll-loops -fomit-frame-pointer  -ffast-math -ftree-vectorize
CXXFLAGS  +=  -g -no-pie -std=c++11 $(COPT)     -Wall -Wextra  -Wno-write-strings  -Wno-deprecated -fopenmp  $(hdrdir)
LDFLAGS +=  $(CXXFLAGS) $(libdir) -lpng -ltiff -ljpeg  -lfftw3 -lfftw3f -fopenmp -lpthread -lX11

LIBMX= ../library/*.o

OBJ	= src_histogram.o src_histogram_midway.o src_histogram_specification.o src_histogram_specification_local.o src_info.o src_linear_combination.o src_gray.o src_gray_4.o src_lp_dist.o src_diff.o src_thre.o src_const.o src_abs.o src_std.o src_mult.o src_mean.o src_get_channels.o src_put_channels.o src_add.o src_convert.o src_extract.o src_yuv_2_rgb.o src_rgb_2_yuv.o src_rotate.o src_float_2_char.o  src_gradient.o src_binarize.o src_iioconvert.o src_put_channels_4.o src_linear_combination_mask.o src_gamma_correction.o src_median.o src_sqrt.o src_conrat_2_rgb.o src_rgb_2_lenzcarmona.o src_rgb_2_conrat.o src_lenzcarmona_2_rgb.o src_psnr.o src_padding.o src_mirror.o src_add_gaussian_noise.o src_add_salt_pepper_noise.o src_paste.o src_rgb_2_cielab.o src_cielab_2_rgb.o src_rgb_2_ycbcr.o src_ycbcr_2_rgb.o src_rgb_2_hsv.o src_hsv_2_rgb.o src_rgb_2_hsl.o src_hsl_2_rgb.o src_rgb_2_hsi.o src_hsi_2_rgb.o src_convert_txt2image.o src_info_height.o src_info_width.o src_sam.o src_histogram_midway_sequence.o src_read_binary_image.o src_equalize_meanstd.o src_equalize_meanstd_local.o src_ssim.o src_pca_svd.o src_pca_color_transform.o src_pca_color_transform_back.o

BIN = src_histogram src_histogram_midway src_histogram_specification src_histogram_specification_local src_info src_linear_combination src_gray src_gray_4 src_lp_dist src_diff src_thre src_const src_abs src_std src_mult src_mean src_get_channels src_put_channels src_add src_convert src_extract src_rgb_2_yuv src_yuv_2_rgb src_rotate src_float_2_char src_gradient src_binarize src_iioconvert src_put_channels_4 src_linear_combination_mask src_gamma_correction src_median src_sqrt src_conrat_2_rgb src_rgb_2_lenzcarmona src_rgb_2_conrat src_lenzcarmona_2_rgb src_psnr src_padding src_mirror src_add_gaussian_noise src_add_salt_pepper_noise src_paste src_rgb_2_cielab src_cielab_2_rgb src_rgb_2_ycbcr src_ycbcr_2_rgb src_rgb_2_hsv src_hsv_2_rgb src_rgb_2_hsl src_hsl_2_rgb src_rgb_2_hsi src_hsi_2_rgb src_convert_txt2image src_info_height src_info_width src_sam src_histogram_midway_sequence src_read_binary_image src_equalize_meanstd src_equalize_meanstd_local src_ssim src_pca_svd src_pca_color_transform src_pca_color_transform_back

LIBBIN=../bin

all: $(OBJ)	$(BIN) 


$(OBJ) : %.o : %.cpp
	$(CXX) -c    $< -o $@ $(CXXFLAGS)


$(BIN) : % : %.o  $(LIBMX)
	$(CXX)   -o $(LIBBIN)/$@  $^ $(LDFLAGS)



.PHONY : clean
clean:
	$(RM) $(OBJ) ; cd $(LIBBIN); rm $(BIN)





