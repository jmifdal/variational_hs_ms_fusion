FROM ubuntu:18.04
WORKDIR /home
COPY . /home/
RUN apt-get update && apt-get install -y make gcc g++\
     libpng-dev libtiff-dev \
     libjpeg-dev libfftw3-dev libx11-dev \
     bc\
     nano 

RUN mkdir bin

RUN cd library && make 
RUN cd libFilter && make 
RUN cd libDomain && make
RUN cd hyperspectral && make
RUN cd libBasic && make

RUN apt-get update

ENV PATH "$PATH:/home/bin"



#COPY --from=build /home/data/fused* $PWD/Desktop/
#RUN echo "export PATH=/home/bin:$PATH"  >> /root/.bashrc


#[FIXME]commande dans le terminal, à enlever plus tard
#for building the image
#docker build -t hsmsfusion .

#for accessing the container with a shell and mounting the volume
#docker run --name demofusion -ti -v $PWD/data:/data hsmsfusion

#generate the hyperspectral and multispectral data with a SNR of 35
#sh batch_hyperspectral_data.sh imgd3 Nikon_D700 35

#launch the fusion. The number "9" corresponds to the weight option used in our algorithm
#sh batch_hyperspectral_fusion.sh imgd3 Nikon_D700 Nikon_D700_transposed 9

#copying the fusion on the host
#Detach from the container without stopping it by pressing CTRL+P followed by CTRL+Q
#recover the ID of the container "demofusion" after typing "docker ps -a"
#then type (in host terminal): sudo docker cp demofusion:/home/demo/fused.tif .

#At this point you have the fusion result in the same folder where you have the docker file
#you can change the location to another one if you wish by replacing "." by the path on your computer 



#for mounting the folder with the data on the container
#docker run -ti -v $PWD/data:/data hsmsfusion





