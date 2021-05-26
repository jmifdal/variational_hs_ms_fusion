##
## First, we compute the NL weights and store all required data as tif files.
##
## Options available:
##
##   - 1 : weights computed on MS (the same for all bands).
##   - 2 : weights computed on MS with Smatrix coefficients.
##   - 3 : weights computed on each upsampled HS band.
##   - 4 : weights computed on full upsampled HS (the same for all bands).
##   - 5 : weights computed on nearest upsampled HS bands.
##   - 6 : neighbors selected on full upsampled HS and weights computed just on associated HS band.
##   - 7 : neighbors selected on full MS and weights computed just on associated HS band.
##   - 8 : neighbors selected on full MS with Smatrix coefficients and weights computed just on associated HS band.
##   - 9 : weights computed on MS with Smatrix coefficients (NEW CODE of the weight option 2 ).
##
## src_hyperspectral_fusion_L1_weights -c $hPatch -d $hDist -n $numNeigh -b $bandSize -r $resSize -p $patchSize $name.multi.tif
##                                     $name.hyper.int.tif $S.tif $weightOpt wxy.tif wyx.tif posxy.tif posyx.tif posw.tif
##                                     numxy.tif numyx.tif
##
## hPatch     filtering parameter for patch distance in NL weights
## hDist	  filtering parameter for spatial distance in NL weights
## numNeigh	  number of (most similar) neighboring pixels used in NL weights
## bandsize   half-size of spectral window in NL weights
## resSize	  half-size of research window in NL weights
## patchSize  half-size of patches in NL weights
##
## name.multi.noisy.tif      input MS noisy image
## name.hyper.noisy.int.tif  input bicubicly interpolated HS noisy image
## S.tif                     input spectral-downsampling matrix
## weightOpt                 input NL weight option (1-9)
## wxy.tif                   output weights wxy
## wyx.tif                   output weights wxy
## posxy.tif                 output neighbour pixel positions posxy
## posyx.tif                 output neighbour pixel positions posyx
## posw.tif                  output weight positions
## numxy.tif                 output number of neighbours per pixel in wxy weights
## numyx.tif                 output number of neighbours per pixel in wyx weights
##
##
## Then, we apply primal-dual fusion algorithm
##
## src_hyperspectral_fusion_L1 -m $lmbM -h $lmbH -u $lmbP -e $tol -i $maxIter -t $primalStep -s $dualStep $name.multi.noisy.tif
##                             $name.hyper.noisy.tif $name.hyper.noisy.int.tif $name.pan.noisy.tif $name.pan.noisy.int.tif wxy.tif wyx.tif
##                             posxy.tif posyx.tif posw.tif numxy.tif numyx.tif $S.tif $St.tif $blur $sampling fused.tif
##
## lmbM	       trade-off parameter for multispectral data term
## lmbH	       trade-off parameter for hyperspectral data term
## lmbP	       trade-off paramater for radiometric constraint
## tol	       stopping precision of primal-dual algorithm
## maxIter	   max allowed iterations for primal-dual algorithm
## primalStep  step-size of primal-dual gradient descent algorithm
## dualStep    step-size of primal-dual gradient ascent algorithm
##
## name.multi.noisy.tif      input MS noisy image
## name.hyper.noisy.tif      input HS noisy image
## name.hyper.noisy.int.tif  input bicubicly interpolated HS noisy image
## name.pan.noisy.tif	     input panchromatic noisy image
## name.pan.noisy.int.tif	 input bicubicly interpolated panchromatic noisy image
## wxy.tif                   input weights wxy
## wyx.tif                   input weights wxy
## posxy.tif                 input neighbour pixel positions posxy
## posyx.tif                 input neighbour pixel positions posyx
## posw.tif                  input weight positions
## numxy.tif                 input number of neighbours per pixel in wxy weights
## numyx.tif                 input number of neighbours per pixel in wyx weights
## S.tif                     input spectral-downsampling matrix
## St.tif                    input transposed spectral-downsampling matrix
## blur                      input standard deviation of Gaussian kernel for spatial downsampling
## sampling                  input sampling factor for spatial downsampling
## fused.tif                 output fused image
##


### Inputs
name=$1        # name of the image (without .tif)
S=$2           # name of the spectral response (without .tif)
St=$3          # name of the transposed spectral response (without .tif)
weightOpt=$4   # choose weight option 1-9


## txt files
fileRMSE=rmse_all.txt
fileRMSEmin=rmse_min.txt
fileSAM=sam_all.txt
fileSAMmin=sam_min.txt
filePSNR=psnr_all.txt
filePSNRmax=psnr_max.txt

### Spatial downsampling parameters
blur=2
sampling=4


## Primal-dual algorithm parameters
tol=0.000001 #the minimum relative error we ask for (the difference between two iterations)
primalStep=0.05 #original step
dualStep=0.05	 #original step

#primalStep=0.0005
#dualStep=0.0005

#maxIter=2000
maxIter=1500


## NL weight parameters
bandSize=1
resSize=7
#patchSize=3
patchSize=1
numNeigh=15
#numNeigh=6
#numNeigh=1;
hDist=2.5

#Limit for sam value below which the result is suspicious
LimitSamError=0.001

## Optimization
rmseMin=$( src_lp_dist -b 10 $name.tif $name.hyper.noisy.int.tif )
samMin=$( src_sam -b 10 $name.tif $name.hyper.noisy.int.tif )
psnrMax=$( src_psnr -b 10 $name.tif $name.hyper.noisy.int.tif )

echo "initialization: psnr = $psnrMax, rmse = $rmseMin, sam = $samMin"
echo "$rmseMin" >> $fileRMSEmin
echo "$samMin" >> $fileSAMmin
echo "$psnrMax" >> $filePSNRmax

hPatch=0.01

	echo "The experiments are running for hPatch = $hPatch"

    src_hyperspectral_fusion_L2_weights -c $hPatch -d $hDist -n $numNeigh -b $bandSize -r $resSize -p $patchSize $name.multi.noisy.tif $name.hyper.noisy.int.tif $S.tif $weightOpt wxy.tif wyx.tif posxy.tif posyx.tif posw.tif numxy.tif numyx.tif
    
    lmbP=0.001
    lmbM=0.1
    lmbH=50
 
    src_hyperspectral_fusion_L2 -m $lmbM -h $lmbH -u $lmbP -e $tol -i $maxIter -t $primalStep -s $dualStep $name.multi.noisy.tif $name.hyper.noisy.tif $name.hyper.noisy.int.tif $name.pan.noisy.tif $name.pan.noisy.int.tif wxy.tif wyx.tif posxy.tif posyx.tif posw.tif numxy.tif numyx.tif $S.tif $St.tif $blur $sampling fused_wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif

    
    rmse=$( src_lp_dist -b 10 $name.tif fused_wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif )
    sam=$( src_sam  -b 10 $name.tif fused_wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif )
    psnr=$( src_psnr  -b 10 $name.tif fused_wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif )

    echo "psnr = $psnr, rmse = $rmse, sam = $sam - hD: $hDist, hP: $hPatch, lP: $lmbP, lM: $lmbM, lH: $lmbH"
    echo "$rmse - hDist: $hDist, hPatch: $hPatch, lmbP: $lmbP, lmbM: $lmbM, lmbH: $lmbH, psnr=$psnr and sam=$sam" >> $fileRMSE
    echo "$sam -  hDist: $hDist, hPatch: $hPatch, lmbP: $lmbP, lmbM: $lmbM, lmbH: $lmbH, rmse=$rmse and psnr=$psnr" >> $fileSAM
    echo "$psnr -  hDist: $hDist, hPatch: $hPatch, lmbP: $lmbP, lmbM: $lmbM, lmbH: $lmbH, rmse=$rmse and sam=$sam" >> $filePSNR

    isSmallerRMSE=$(echo "$rmse<=$rmseMin" | bc)
    isSmallerSAM=$(echo "$sam<=$samMin" | bc)
    isGreaterPSNR=$(echo "$psnr>=$psnrMax" | bc)
    isGreaterSAMthanLimit=$(echo "$sam>=$LimitSamError" | bc)


    if [ "$isSmallerRMSE" -eq 1 ]; then
        echo "$rmse - hP: $hPatch, lP: $lmbP, lM: $lmbM, lH: $lmbH, psnr=$psnr and sam=$sam" >> $fileRMSEmin
        #cp fused.tif fused.wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif
        rmseMin=$rmse
    fi

    if [ "$isSmallerSAM" -eq 1 ] && [ "$isGreaterSAMthanLimit" -eq 1 ]; then
        echo "$sam - hP: $hPatch, lP: $lmbP, lM: $lmbM, lH: $lmbH, rmse=$rmse and psnr=$psnr" >> $fileSAMmin
        #cp fused.tif fused.wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif
         samMin=$sam
    fi

    if [ "$isGreaterPSNR" -eq 1 ]; then
        echo "$psnr - hP: $hPatch, lP: $lmbP, lM: $lmbM, lH: $lmbH, rmse=$rmse and sam=$sam" >> $filePSNRmax
        #cp fused.tif fused.wOpt$weightOpt.hP$hPatch.lP$lmbP.lM$lmbM.lH$lmbH.tif
        psnrMax=$psnr
    fi




 exit
