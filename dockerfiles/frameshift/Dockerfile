FROM ubuntu:latest

MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)
ENV DEBIAN_FRONTEND noninteractive


RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        python3.8-minimal \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

WORKDIR /framecorrect/
COPY correct_frameshift.py /framecorrect/

ENV LANG pt_BR.UTF-8