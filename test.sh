#!/bin/bash

for i in `seq 1 10`;
do
    ./run.sh $1 $2 $3 results/$4
done