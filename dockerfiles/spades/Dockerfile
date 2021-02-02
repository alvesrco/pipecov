FROM ubuntu:groovy
MAINTAINER ALVES Ronnie <alvesrco@gmail.com>
#Updated by Renato Oliveira 

RUN apt-get -y update
RUN apt-get -y dist-upgrade 
RUN apt-get install --no-install-recommends --no-install-suggests -y python wget

#RUN ln -s /usr/bin/python3 /usr/bin/python 

	
#RUN pip3 install regex

#RUN  wget http://cab.spbu.ru/files/release3.15.0/SPAdes-3.15.0-Linux.tar.gz -O /tmp/spades.tar.gz && \
#     tar -xvzf /tmp/spades.tar.gz -C /opt/ && \
#     rm /tmp/spades.tar.gz

RUN  wget http://cab.spbu.ru/files/release3.15.0/SPAdes-3.15.0.tar.gz -O /tmp/spades.tar.gz && \
     tar -xvzf /tmp/spades.tar.gz -C /opt/ && \
     rm /tmp/spades.tar.gz

RUN apt-get install --no-install-recommends --no-install-suggests -y cmake make g++ zlib1g-dev libbz2-dev bzip2 libz-dev

RUN apt-get -y autoremove
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

RUN cd /opt/SPAdes-3.15.0/ && ./spades_compile.sh

ENV PATH="/opt/SPAdes-3.15.0/:${PATH}"

WORKDIR /data/
