# Steps to clone the baydem github repository and run the baydem tests inside
# a Docker container. Really, the only reason the first clone command is needed
# is to obtain this Dockerfile and setup.R. When running the Docker container,
# a mirrored directory is used, which is not necessary for just running the
# tests, but is usually necessary for running analyses. Modify the Docker
# commands to fit your situation.
#
# git clone https://github.com/eehh-stanford/baydem
# cd baydem
# docker build -t michaelholtonprice/baydem .
# docker run --name baydem -itv //c/mirrored_baydem_data:/data michaelholtonprice/baydem
# git clone https://github.com/eehh-stanford/baydem
# cd baydem
# R
# library(devtools)
# test()

FROM ubuntu:20.04

# Set the following environmental variable to avoid interactively setting the
# timezone with tzdata when installing R
ENV DEBIAN_FRONTEND=noninteractive

# Install R version 4
RUN apt-get update && \
    apt-get install -y gnupg2 &&\
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    apt-get install -y software-properties-common && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update && \
    apt-get install -y r-base &&\
    apt-get install -y vim && \
    apt-get install -y git && \
#    apt-get install apt-transport-https && \
#    apt-get install software-properties-common && \
    apt-get clean

# Install the tools needed for subsequent R package installations
RUN apt-get install -y libssl-dev && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libfontconfig1-dev && \
    apt-get install -y libharfbuzz-dev && \
    apt-get install -y libfribidi-dev && \
    apt-get install -y libfreetype6-dev && \
    apt-get install -y libpng-dev && \
    apt-get install -y libtiff5-dev && \
    apt-get install -y libjpeg-dev && \
    apt-get install -y libv8-dev

# The following line can be copied at the command line to install the preceding
# requirements:
# sudo apt-get install libssl-dev libxml2-dev libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libv8-dev

COPY setup.R .
RUN Rscript setup.R

