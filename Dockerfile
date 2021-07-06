FROM coffeateam/coffea-base:latest

WORKDIR /home/cmsusr/

RUN pip install dask-jobqueue distributed
