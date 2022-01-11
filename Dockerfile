FROM ubuntu:bionic
LABEL org.opencontainers.image.authors=nicolasbock@gmail.com

COPY scripts/prepare-container.sh /usr/sbin
RUN /usr/sbin/prepare-container.sh
COPY scripts/install-bml.sh /usr/sbin

ENV CC gcc-10
ENV CXX g++-10
ENV FC gfortran-10
ENV CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Debug}
RUN INSTALL_DIR=/usr /usr/sbin/install-bml.sh

WORKDIR /root
