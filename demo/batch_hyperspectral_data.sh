##
## Generate data needed for multispectral and hyperspectral fusion.
##
## Gaussian noise is added to the generated multispectral and hyperspectral data.
##
## - $name.multi.tif           : high-resolution MS image.
## - $name.multi.noisy.tif     : high-resolution MS noisy image.
## - $name.hyper.tif           : low-resolution HS image.
## - $name.hyper.noisy.tif     : low-resolution HS noisy image.
## - $name.hyper.noisy.int.tif : bicubicly interpolated HS noisy image.
## - $name.pan.noisy.tif       : panchromatic noisy image obtained by averaging MS noisy channels.
## - $name.pan.noisy.int.tif   : bicubicly interpolated panchromatic noisy image.
##


### Inputs
name=$1    # name of the image (without .tif)
S=$2       # name of the spectral response (without .tif)
snr=$3     # signal to noise ratio of the noise


### Spatial downsampling parameters
blur=2
sampling=4


## Multispectral image
src_hyperspectral_generate_multispectral $name.tif $name.multi.tif $S.tif
src_add_gaussian_noise -r $snr $name.multi.tif $name.multi.noisy.tif


## Hyperspectral image
src_filter_convol_gaussian -b 1 $name.tif tmp.tif $blur
src_sample tmp.tif $name.hyper.tif $sampling
src_add_gaussian_noise -r $snr $name.hyper.tif $name.hyper.noisy.tif


## Interpolated hyperspectral image
src_spline_zoom -z $sampling -v $name.hyper.noisy.tif $name.hyper.noisy.int.tif


## Panchromatic image
src_hyperspectral_generate_hyperpan $name.multi.noisy.tif $name.pan.noisy.tif $S.tif


## Interpolated panchromatic image
src_filter_convol_gaussian -b 1 $name.multi.noisy.tif tmp.tif $blur
src_sample tmp.tif tmp.tif $sampling
src_spline_zoom -z $sampling -v tmp.tif tmp.tif
src_hyperspectral_generate_hyperpan tmp.tif $name.pan.noisy.int.tif $S.tif


rm tmp.tif
