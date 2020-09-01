FROM ubuntu:16.04

MAINTAINER Roberto Xavier <rbxjunior@gmail.com>

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN apt-get install -y wget vim dpkg
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.471-1_amd64.deb
RUN dpkg -i mafft_7.471-1_amd64.deb

WORKDIR /data
