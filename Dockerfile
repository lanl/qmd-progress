FROM ubuntu:bionic

COPY prepare-container.sh /usr/sbin
RUN /usr/sbin/prepare-container.sh
COPY scripts/install-bml.sh /usr/sbin
RUN /usr/sbin/install-bml.sh

WORkDIR /root
