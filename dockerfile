FROM ubuntu

RUN apt update -y && apt upgrade -y
RUN apt install python3 python3-pip -y 
RUN apt-get install -y dotnet6
RUN apt-get install -y dotnet-runtime-6.0

WORKDIR /P6

COPY ./ ./

RUN python3 setup.py


