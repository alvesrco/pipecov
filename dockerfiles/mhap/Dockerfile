FROM openjdk:8-jdk-alpine

LABEL Roberto Xavier <rbxjunior@gmail.com>

RUN apk add --no-cache curl tar bash procps wget

# Install Maven
RUN wget http://us.mirrors.quenda.co/apache/maven/maven-3/3.6.3/binaries/apache-                                                                                                                     maven-3.6.3-bin.tar.gz -P /tmp
RUN tar xzvf /tmp/apache-maven-3.6.3-bin.tar.gz -C /root
ENV PATH /root/apache-maven-3.6.3/bin:$PATH

# Install MHAP
RUN wget https://github.com/marbl/MHAP/archive/v2.1.1.tar.gz -P /tmp
RUN tar -xzf /tmp/v2.1.1.tar.gz -C /root
RUN cd /root/MHAP-2.1.1 && mvn install

WORKDIR /root/MHAP-2.1.1/target
