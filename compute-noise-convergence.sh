#!/bin/bash

nb=$#
args=$*
echo "#args=${nb}"
echo " args=${args}"
if test ${nb} -lt 7; then 
    echo "Usage: $0 estimator shape kernel R r alpha noise"
    exit 0
fi
estimator=$1
shape=$2
kernel=$3
bigR=$4
smallR=$5
alpha=$6
noiselvl=$7
h_begin=0.75
h_mult=0.5
h_end=0.04
min_bbox=-10
max_bbox=10
trivial_t=8

case ${estimator} in
"VCM")     ;;
*)         echo "Bad estimator: ${estimator}"; exit 1;;
esac

case ${shape} in
"ellipse") polynomial="90-3*x^2-2*y^2-z^2";;
"rcube")   polynomial="6561-x^4-y^4-z^4";;
"goursat") polynomial="8-0.03*x^4-0.03*y^4-0.03*z^4+2*x^2+2*y^2+2*z^2";;
"distel")  polynomial="10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))";;
"leopold") polynomial="100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)";;
"diabolo") polynomial="x^2-(y^2+z^2)^2";;
"heart")   polynomial="-1*(x^2+2.25*y^2+z^2-1)^3+x^2*z^3+0.1125*y^2*z^3";;
"crixxi")  polynomial="-0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3";;
"torus")   polynomial="-1*(x^2+y^2+z^2+6*6-2*2)^2+4*6*6*(x^2+y^2)";;
*)         echo "Bad shape: ${shape}"; exit 1;;
esac

case ${kernel} in
"hat")     ;;
"ball")     ;;
*)         echo "Bad kernel: ${kernel}"; exit 1;;
esac


h=${h_begin}
while test `echo "${h} >= ${h_end}"|bc -l` -eq 1; do
    file="${estimator}"
    if test ! -d ${file}; then mkdir ${file}; fi
    file="${estimator}/noise-${noiselvl}"
    if test ! -d ${file}; then mkdir ${file}; fi
    file="${estimator}/noise-${noiselvl}/${kernel}-${alpha}"
    if test ! -d ${file}; then mkdir ${file}; fi
    file="${estimator}/noise-${noiselvl}/${kernel}-${alpha}/${shape}"
    toExport=""
    if test `echo "(${max_bbox}-(${min_bbox}))/${h} < 512"| bc -l` -eq 1; then \
        toExport="-x ${estimator}"; \
    fi
    #toExport="-x VCM"
    echo "------------------------------------------------------------------"
    echo "#      file      =${file}"
    echo "#      estimator =${estimator}"
    echo "#      shape     =${shape}"
    echo "#      polynomial=${polynomial}"
    echo "#      kernel    =${kernel}"
    echo "#      gridstep  =${h}"
    echo "#      noiselvl  =${noiselvl}"
    echo "#"
    echo "./implicitShape3NormalEstimation -p ${polynomial} -o ${file} -a ${min_bbox} -A ${max_bbox} -e ${estimator} -R ${bigR} -r ${smallR} -t ${trivial_t} -E 1 ${toExport} --alpha ${alpha} -k ${kernel} -N ${noiselvl} -g ${h}"
    ./noisyImplicitShape3NormalEstimation -p ${polynomial} -o ${file} -a ${min_bbox} -A ${max_bbox} -e ${estimator} -R ${bigR} -r ${smallR} -t ${trivial_t} -E 1 ${toExport} --alpha ${alpha} -k ${kernel} -N ${noiselvl} -g ${h}
    echo "------------------------------------------------------------------"
    echo
    echo
    h=0`echo "${h}*${h_mult}"|bc -l`
done
