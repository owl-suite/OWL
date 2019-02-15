# Script for downloading and compiling Quantum Espresso from scratch


#!/bin/bash

  # Download espresso-5.4.0 if it does not exist
  if [ ! -d q-e-qe-5.4 ]; then
    wget https://github.com/QEF/q-e/archive/qe-5.4.tar.gz
    tar zxvf qe-5.4
  fi

  # Build CPU code
  #if [ "$1" == "cpu" ]; then
    cd q-e-qe-5.4
    ./configure --enable-openmp CC=mpicc F77=mpif77 MPIF90=mpif90
    make pw
    #make all
    cd ..
    rm ../qe
    ln -s q-e-qe-5.4 ../qe
  #fi

# Build GPU code
#  if [ "$1" == "gpu" ]; then
#    # Download QE-GPU-5.4.0.tar.gz if it does not exist
#    #if [ ! -f QE-GPU-5.4.0.tar.gz ]; then 
#    #  wget http://www.qe-forge.org/gf/download/frsrelease/211/972/QE-GPU-5.4.0.tar.gz; 
#    #fi
#   
#    cp QE-GPU-5.4.0.tar.gz q-e-qe-5.4
#    mv q-e-qe-5.4 q-e-qe-5.4_gpu
#    cd q-e-qe-5.4_gpu
#    tar zxvf QE-GPU-5.4.0.tar.gz
#    cd GPU
#    ./configure --enable-parallel --enable-openmp --enable-cuda \
#    --with-gpu-arch=sm_35 --with-cuda-dir=$CRAY_CUDATOOLKIT_DIR \
#    --with-phigemm  --with-scalapack --with-elpa \
#    --without-magma CC=mpicc F77=mpif77 MPIF90=mpif90
#    cd ..
#    make -f Makefile.gpu pw-gpu
#  fi

