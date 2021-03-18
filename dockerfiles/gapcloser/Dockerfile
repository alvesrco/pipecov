FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y install --no-install-recommends --no-install-suggests wget
RUN apt-get -y install --no-install-recommends --no-install-suggests unzip
RUN apt-get -y autoremove; \
	apt-get clean; \
	rm -rf /var/lib/apt/lists/*


ENV ZIP=GapCloser-bin-v1.12-r6.tgz
ENV URL=https://github.com/BioContainers/software-archive/releases/download/gapcloser
ENV DST=/usr/local/bin

RUN wget --no-check-certificate $URL/$ZIP -O $DST/$ZIP
RUN tar xvf $DST/$ZIP -C $DST
RUN rm $DST/$ZIP

COPY gapcloser.config /data/

WORKDIR /data/