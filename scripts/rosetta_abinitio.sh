#!/bin/bash

main_dir=~/psp

protein=$1
it=$2

mkdir -p $main_dir/results/$protein/rosetta
cd $main_dir/results/$protein/rosetta

mkdir -p tmp_$it
cd tmp_$it

$main_dir/external/rosetta/bin/AbinitioRelax.linuxgccrelease -in:file:fasta $main_dir/proteins/$protein/fasta -in:file:frag3 $main_dir/proteins/$protein/frag3 -in:file:frag9 $main_dir/proteins/$protein/frag9 -out:pdb -abinitio::increase_cycles 28

cp *.pdb ../decoy_$it.pdb

cd ..
rm -rf tmp_$it
