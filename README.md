# baydem
The R package baydem (for bayESIAN demOGRAPHY) provides tools for reconstructing past and present demography.

# Installation
There are three options for installing and using the package:

(1) Install on an existing computer with necessary dependencies
(2) Build a Docker image using the Dockerfile
(3) Use a Docker image from Docker Hub

## Option 1: Install on an existing computer with necessary dependencies
Details will vary based on the machine and operating system since certain R packages on which baydem depends require the installation of tools outside R. The following should be sufficient:

(a) Install devtools with its dependencies
https://www.rdocumentation.org/packages/devtools

(b) Install rstan with its dependencies
http://mc-stan.org/rstan/

(c) Install baydem
```
library(devtools)
install_github("eehh-stanford/baydem")
```

## Option 2: Build a Docker image using the Dockerfile
First, clone the database and build the Docker image.
```bash
git clone https://github.com/eehh-stanford/baydem
cd baydem
docker build -t michaelholtonprice/baydem .
```

Start a container. Use the -v tag to mirror a directory for passing files between the host machine and the Docker container. The directory to the left of the semicolon is for the host machine and the directory to the right of the semicolon is for the Docker container. The path for the host machine will need to be modified for your situation.

```bash
docker run --name baydem -itv //c/mirrored_baydem_data:/data michaelholtonprice/baydem
```

## Option 3: Use a Docker image from Docker Hub
TODO: update example