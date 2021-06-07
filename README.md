# Variational fusion of hyperspectral and multispectral images


## Docker

This project is coded with C++ and the libraries needed for its functionning are included in a Dockerfile which enables our code to work on any OS based systems

* First Docker needs to be installed on you system with the help of this link: https://docs.docker.com/engine/install/

* Clone our project in the location of your choosing on your local system and open a terminal in the cloned folder

## Building the docker container

From the Dockerfile we're going to build an image called "hsmsfusion" in the folder downloaded from our Github repository

```bash
docker build -t hsmsfusion .
```

## Fusion

At this point all the libraries and the packages are installed. The next commannd line does multiple manouvers: it starts a docker container named ```demofusion``` from the built image, it mounts a volume on the docker container that points to the ```demo``` folder in downloaded one and launches the ```demo.sh``` script. 

The ```demo.sh``` script executes three commands: it moves to the ```demo``` folder, it creates the HS and MS images and all the data needed for the fusion and finally, it carries out the fusion with the generated data. 

```bash
docker run --name demofusion -v $PWD/demo:/home/demo hsmsfusion sh demo.sh
```
## Recovering the fusion result
The fusion result will be available in the ```demo``` folder in the cloned repository.  
