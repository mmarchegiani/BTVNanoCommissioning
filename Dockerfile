FROM coffeateam/coffea-dask-cc7:latest
#FROM coffeateam/coffea-base:latest

WORKDIR /home/cmsusr/

COPY . /home/cmsusr/

#RUN pip install dask-jobqueue distributed psutil
#RUN apt update && \
#    apt upgrade -y && \
#    apt install strace -y
#RUN yum install strace -y 
