FROM ubuntu

RUN apt update && apt upgrade && apt install python3 python3-pip -y 
WORKDIR /P6

COPY . .
RUN pip install -r requirements.txt