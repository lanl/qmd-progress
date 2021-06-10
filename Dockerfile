FROM ubuntu:bionic

COPY prepare-container.sh /usr/sbin
RUN /usr/sbin/prepare-container.sh
COPY scripts/install-bml.sh /usr/sbin

ENV CC gcc-6
ENV CXX g++-6
ENV FC gfortran-6

RUN /usr/sbin/install-bml.sh

WORkDIR /root
