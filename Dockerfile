FROM ubuntu:latest

RUN apt-get update && \
    apt-get install gfortran -y

WORKDIR /app
