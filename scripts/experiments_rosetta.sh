#!/bin/bash

proteins=(1vii  1gb1  1hhp  1i6c  1crn  2imu  T0868  T0900)

for protein in ${proteins[@]};
do
    args=""
    for i in `seq 1 32`;
    do
	    args+="$protein $i "
    done

    echo $args | xargs -n 2 -P 32 ./rosetta_abinitio.sh
done
