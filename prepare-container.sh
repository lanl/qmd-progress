#!/bin/bash

set -e -u -x

SUDO=$(which sudo || true)

for i in $(seq 5); do
  ${SUDO} apt-get update && break
done

${SUDO} apt-get install --assume-yes --no-install-recommends \
  apt-transport-https \
  ca-certificates \
  gnupg \
  wget

cat <<EOF | ${SUDO} tee /etc/apt/sources.list.d/toolchain.list
deb http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu bionic main
# deb-src http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu bionic main
EOF
${SUDO} apt-key adv --keyserver keyserver.ubuntu.com \
  --recv-keys 60C317803A41BA51845E371A1E9377A2BA9EF27F

for i in $(seq 5); do
  ${SUDO} apt-get update && break
done

${SUDO} ln -fs /usr/share/zoneinfo/UTC /etc/localtime
${SUDO} apt-get install --assume-yes tzdata
DEBIAN_FRONTEND=noninteractive ${SUDO} dpkg-reconfigure \
  --frontend noninteractive tzdata

${SUDO} apt-get install --assume-yes --no-install-recommends \
  cmake \
  cmake-data \
  g++-6 \
  gcc-6 \
  gfortran-6 \
  libblas-dev \
  liblapack-dev \
  libmetis-dev \
  libopenmpi-dev \
  make \
  pkg-config \
  python \
  python3-numpy \
  python-numpy
