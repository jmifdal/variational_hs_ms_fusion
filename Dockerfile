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