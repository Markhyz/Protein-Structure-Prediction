#!/bin/bash

#proteins=(1vii  1gb1  1hhp  1i6c  1crn  2imu  T0868  T0900)
#proteins=(1acw  1ail  1crn  1enh  1rop  1zdd  2mr9  2p81)
proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

for protein in ${proteins[@]};
do
    args=""
    for i in `seq 1 32`;
    do
	    args+="$protein $i "
    done

    echo $args | xargs -n 2 -P 32 ./scripts/rosetta_abinitio.sh
done
