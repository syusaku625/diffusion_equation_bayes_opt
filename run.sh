#!/bin/bash

EXE=/home/osin1112/diffusion_equation_bayes_opt/out/Main
DIR=/home/osin1112/diffusion_equation_bayes_opt/example/
TP=${DIR}test.tp
NO=11

cd ${DIR}No.${NO}/normal
${EXE} ${TP} > opt.log

#OUT="tmp"
#while [ -n ${OUT} ]
#do
#    echo loop
#    rm opt.log
#    OUT=`${EXE} ${TP} > opt.log`
#done

#for i in 21 22 23 24 25 26 27 28 29 30
#do
#    cd ${DIR}No.${NO}/phi_vessel_max_${i}
#    ${EXE} ${TP} > opt.log
#done