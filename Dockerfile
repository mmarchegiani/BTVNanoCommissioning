FROM coffeateam/coffea-base:latest

WORKDIR /home/cmsusr/

# Add code
COPY . /home/cmsusr/
#RUN git clone https://gitlab.cern.ch/algomez/BTVNanoCommissioning.git -b v01
