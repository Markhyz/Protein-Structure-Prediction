#!/bin/bash

proteins=(1acw  1ail  1crn  1enh  1i6c  1rop  1zdd  2mr9  2p81  T0868  T0900)

for protein in ${proteins[@]};
do
	mkdir -p ./statistics/$1/$protein
	cp ./results/$protein/best_ca_rmsd2 ./statistics/$1/$protein/rmsd
	cp ./results/$protein/best_ca_gdt2 ./statistics/$1/$protein/gdt
	cp ./results/$protein/mufold_rmsd ./statistics/$1/$protein/mufold_rmsd
	cp ./results/$protein/mufold_gdt ./statistics/$1/$protein/mufold_gdt
	cp ./results/$protein/rosetta_rmsd ./statistics/$1/$protein/rosetta_rmsd
	cp ./results/$protein/rosetta_gdt ./statistics/$1/$protein/rosetta_gdt
done

