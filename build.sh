#!/bin/sh
mkdir build
cd build
cmake -D bayesopt_DIR=/mnt/c/lib/bayesopt \
-D TP_DIR=/mnt/c/share/lib/TextParser \
-D CMAKE_INSTALL_PREFIX=/home/maeda/openga_example \
-D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
..

make -j 4 && make install
