FROM python:2.7.18

# Donwnload and install RNA
RUN wget https://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz
RUN tar -zxvf RNAstructureLinuxTextInterfaces64bit.tgz
RUN cd RNAstructure && \
    make all

ENV PATH=$PATH:/RNAstructure/exe
ENV DATAPATH=/RNAstructure/data_tables

# Setup Metaplotlib
RUN python -m pip install -U pip setuptools
RUN python -m pip install matplotlib

# Setup httplib2
RUN pip install httplib2

# Copy required files
RUN mkdir /Superfold
WORKDIR Superfold
COPY . /Superfold/
