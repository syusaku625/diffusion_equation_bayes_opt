## 三相脳循環モデル+Bayes optimization ![ソースコードサイズ](https://img.shields.io/github/repo-size/syusaku625/diffusion_equation_bayes_opt)

## 依存ライブラリ
![ソースコードサイズ](https://img.shields.io/badge/cmake-v.3.16.3-yellowgreen)
![ソースコードサイズ](https://img.shields.io/badge/boost-v.1.71.0-green)
![ソースコードサイズ](https://img.shields.io/badge/HDF5-v.1.12.2-orange)
![ソースコードサイズ](https://img.shields.io/badge/TextParser-up__to__date-red)
![ソースコードサイズ](https://img.shields.io/badge/bayesopt-up__to__date-blue)
- cmake (https://cmake.org/)
- Boost (https://www.boost.org/)
- HDF5 (https://www.hdfgroup.org/solutions/hdf5/)
- TextParser (https://github.com/avr-aics-riken/TextParser)
- bayesopt (https://github.com/rmcantin/bayesopt)


## How to build
1. please edit build.sh

```
#!/bin/sh
mkdir build
cd build
cmake -D bayesopt_DIR=/path/to/your/bayesopt \
-D TP_DIR=/path/to/your/TextParser \
-D CMAKE_INSTALL_PREFIX=/path/to/your_install_directory \
-D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
..

make -j 4 && make install

```
2. please run build.sh
```
sh build.sh
```

## How to run
please run executable file in example directory  

Give test.tp as command line argument

```
./../Main test.tp
```