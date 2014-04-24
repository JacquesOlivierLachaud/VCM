#!/bin/bash

nb=$#
args=$*
echo "#args=${nb}"
echo " args=${args}"
if test ${nb} -lt 5; then 
    echo "Usage: $0 estimator shape R r alpha"
fi
estimator=$1
shape=$2
bigR=$3
smallR=$4
alpha=$5
h_begin=1
h_mult=0.5
h_end=0.03
min_bbox=-10
max_bbox=10

case ${estimator} in
"VCM") ;;
*)     echo "Bad estimator: ${estimator}";;
esac

case ${shape} in
"ellipse") polynomial="90-3*x^2-2*y^2-z^2";;
"rcube")   polynomial="6561-x^4-y^4-z^4";;
"goursat") polynomial="8-0.03*x^4-0.03*y^4-0.03*z^4+2*x^2+2*y^2+2*z^2";;
*)         echo "Bad shape: ${shape}";;
        exit 1;
esac

./implicitShape3NormalEstimation -p "90-3*x^2-2*y^2-z^2" -o "VCM/ellipse" -a -10 -A 10 -e VCM -R 3 -r 3 -t 2 -E 0 -x "vcm" -g 0.030625
