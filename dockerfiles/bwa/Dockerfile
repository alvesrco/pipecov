FROM alpine:3.11
RUN apk add --no-cache bash gcc musl-dev make bzip2 zlib-dev
COPY bwa-0.7.17.tar.bz2 /tmp/
RUN tar xjf /tmp/bwa-0.7.17.tar.bz2 -C /root
RUN make -C /root/bwa-0.7.17
ENV PATH /root/bwa-0.7.17:$PATH
RUN mkdir /data
WORKDIR /data