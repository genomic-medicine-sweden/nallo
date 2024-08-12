# syntax=docker/dockerfile:1
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y procps sed zstd && rm -rf /var/lib/apt/lists/*
# install app
COPY fqcrs /usr/bin/fqcrs

# final configuration
ENV FLASK_APP=fqcrs
