#!/bin/sh
mkdir build
cd build
cmake -D bayesopt_DIR=/home/osin1112/lib/bayesopt \
-D TP_DIR=/home/osin1112/lib/TextParser \
-D CMAKE_INSTALL_PREFIX=/home/osin1112/diffusion_equation_bayes_opt/out \
-D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
..

make -j 4 && make install
