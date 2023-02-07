#!/bin/bash

set -e -u -x

SUDO=$(which sudo || true)

for i in $(seq 5); do
  ${SUDO} apt-get update && break
done

DEBIAN_FRONTEND=noninteractive ${SUDO} apt-get install --assume-yes --no-install-recommends \
  apt-transport-https \
  ca-certificates \
  gnupg \
  wget

cat <<EOF | ${SUDO} tee /etc/apt/sources.list.d/toolchain.list
deb http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu jammy main
# deb-src http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu jammy main
EOF
${SUDO} apt-key adv --keyserver keyserver.ubuntu.com \
  --recv-keys 60C317803A41BA51845E371A1E9377A2BA9EF27F

cat <<EOF | ${SUDO} tee /etc/apt/sources.list.d/emacs.list
deb http://ppa.launchpad.net/kelleyk/emacs/ubuntu jammy main
# deb-src http://ppa.launchpad.net/kelleyk/emacs/ubuntu jammy main
EOF
${SUDO} apt-key adv --keyserver keyserver.ubuntu.com \
  --recv-keys 873503A090750CDAEB0754D93FF0E01EEAAFC9CD

for i in $(seq 5); do
  ${SUDO} apt-get update && break
done

${SUDO} ln -fs /usr/share/zoneinfo/UTC /etc/localtime
DEBIAN_FRONTEND=noninteractive ${SUDO} apt-get install --assume-yes tzdata
DEBIAN_FRONTEND=noninteractive ${SUDO} dpkg-reconfigure \
  --frontend noninteractive tzdata

DEBIAN_FRONTEND=noninteractive ${SUDO} apt-get install --assume-yes --no-install-recommends \
  build-essential \
  bundler \
  cmake cmake-data \
  emacs \
  gcc-9 g++-9 gfortran-9 \
  gcc-10 g++-10 gfortran-10 \
  gcc-11 g++-11 gfortran-11 \
  git-core \
  indent \
  libblas-dev \
  liblapack-dev \
  libmetis-dev \
  libopenmpi-dev \
  make \
  pkg-config \
  python3 \
  python3-numpy \
  sudo
