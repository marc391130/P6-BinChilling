FROM ubuntu:22.04


RUN apt update -y && apt upgrade -y
RUN apt install python3 python3-pip -y
RUN apt install git -y

RUN apt-get install -y dotnet6
RUN apt-get install -y dotnet-runtime-6.0

WORKDIR /P6

RUN git clone https://github.com/marc391130/P6-BinChilling.git .

RUN python3 setup.py


