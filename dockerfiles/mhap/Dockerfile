FROM openjdk:8-jdk-alpine

LABEL Roberto Xavier <rbxjunior@gmail.com>

RUN apk add --no-cache curl tar bash procps wget

# Install Maven
RUN wget http://us.mirrors.quenda.co/apache/maven/maven-3/3.6.3/binaries/apache-maven-3.6.3-bin.tar.gz -P /tmp
RUN tar xzvf /tmp/apache-maven-3.6.3-bin.tar.gz -C /root
ENV PATH /root/apache-maven-3.6.3/bin:$PATH

# Install MHAP
RUN wget https://github.com/marbl/MHAP/releases/download/v2.1.1/mhap-2.1.1.jar.gz
RUN gunzip mhap-2.1.1.jar.gz

WORKDIR /data
