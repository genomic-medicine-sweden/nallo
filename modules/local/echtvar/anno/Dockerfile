# syntax=docker/dockerfile:1
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y procps sed wget && rm -rf /var/lib/apt/lists/*

# install app
RUN wget -O /usr/bin/echtvar https://github.com/brentp/echtvar/releases/download/v0.2.2/echtvar \
    && chmod +x /usr/bin/echtvar

CMD [ "echtvar" ]
