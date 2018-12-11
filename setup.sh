#!/bin/bash

source /usr/local/sklib_gcc4.8.5/bashenv_gcc4.8.5_skofl_18a+atmpd_18a
export FITQUN_ROOT=${ATMPD_ROOT}/src/recon/fitqun/

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SKOFL_ROOT}/lib:${ATMPD_ROOT}/lib:`root-config --libdir`
