FROM coffeateam/coffea-base:latest

WORKDIR /home/cmsusr/

RUN pip install dask-jobqueue distributed
# Add code
#COPY . /home/cmsusr/
#RUN git clone https://gitlab.cern.ch/algomez/BTVNanoCommissioning.git -b v01
#RUN ls -lR /home/cmsusr/
