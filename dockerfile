FROM ubuntu

RUN apt update && apt upgrade && apt install python3 python3-pip -y 

RUN pip install numpy
RUN pip install tqdm
RUN pip install matplotlib
RUN pip install scipy

WORKDIR /P6

COPY ./ACE ./ACE