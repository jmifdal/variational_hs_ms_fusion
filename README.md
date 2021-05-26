# Variational fusion of hyperspectral and multispectral images


## Docker

This project is coded with C++ and the libraries needed for its functionning are included in a Dockerfile which enables our code to work on any OS based systems

* First Docker needs to be installed on you system with the help of this link: https://docs.docker.com/engine/install/

* Clone our project in the location of your choosing on your local system and open a terminal

## Building the docker container

* From the Dockerfile we're going to build an image called "hsmsfusion" in the folder downloaded from our Github repository

```bash
docker build -t hsmsfusion .
```

* At this point all the libraries and the packages are installed. Now we're going to start a docker container named "demofusion" from the built image and mount a volume on it that points to the downloaded folder

```bash
docker run --name demofusion -ti -v $PWD/demo:/demo hsmsfusion
```

Now we should be in a ubuntu shell and more specifically in the '/home' directory. By typing "ls" we should see the same files as the ones downloaded from the Github

## Fusion
### Data generation

First we generate the HS, MS and all the data that we need for the    fusion with an SNR (Signal to Noise Ratio) of 35. For this we do as follows
```bash
sh batch_hyperspectral_data.sh bicycle Nikon_D700 35
```
### Data fusion
Now we are ready to start fusing!
```bash
sh batch_hyperspectral_fusion.sh bicycle Nikon_D700 Nikon_D700_transposed 9
```
## Recovering the fusion reuslt
Now that our fused image called "fused.tif" is ready, we will bring it back to our local system. 
* First we leave the docker container without stopping it by pressing CTRL+P followed by CTRL+Q. Now we are back to our local terminal.

* Next we copy 'fused.tif' to the same downloaded folder by typing

    ```bash
    sudo docker cp demofusion:/home/demo/fused.tif .
    ```
* At this point we are done and we can stop the contained. We display the container's ID by typing  ```docker ps -a``` and then we stop with ```docker stop ID``` followed by ```docker rm ID```
