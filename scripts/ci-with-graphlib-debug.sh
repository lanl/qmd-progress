#!/bin/bash

set -e -u -x

basedir=$(readlink --canonicalize $(dirname ${BASH_SOURCE[0]})/..)

[[ -f ${basedir}/scripts/ci-defaults.sh ]] && . ${basedir}/scripts/ci-defaults.sh

export PROGRESS_GRAPHLIB=yes
export OMP_NUM_THREADS=4
export COMMAND=testing
