

cd demo

#generation of the HS and MS and imgaes and all the data needed for the fusion
sh batch_hyperspectral_data.sh bookshelves Nikon_D700 35

#fusion
sh batch_hyperspectral_fusion.sh bookshelves Nikon_D700 Nikon_D700_transposed 9

#remove temporary files
rm num* pos* psnr* sam* rmse* wyx.tif wxy.tif 
