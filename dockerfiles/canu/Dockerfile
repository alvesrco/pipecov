FROM ubuntu:16.04

LABEL Roberto Xavier <rbxjunior@gmail.com>

RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install build-essential wget

RUN apt-get install -y curl

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# install canu

RUN wget https://github.com/marbl/canu/archive/v1.8.tar.gz -P /tmp
RUN tar -xzf /tmp/v1.8.tar.gz -C /root
RUN cd /root/canu-1.8/src && make

ENV PATH /root/canu-1.8/Linux-amd64/bin:$PATH

WORKDIR /data
