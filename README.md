# Overview

BRKGA optimizer for the Protein Structure Prediction problem.

# Build

To build this project, you will need to have the source code of [Rosetta](https://www.rosettacommons.org/). The default path to link the source code in this project is `external/rosetta`.

1. Clone this repo
2. Clone the submodules (`git submodule update --init`)
3. Compile the project
    - This project was built using CMake. A CMake file is expected in the Rosetta source directory. You can build a CMake file for Rosetta by either reusing one of their CMake files (I tried without success), or create your own custom Rosetta CMake file.

# Project structure

## Main files

The main program is contained in the `main.cpp` file, which will include the other libraries. The Rosetta libraries are contained in `rosetta.hpp`.

## Evolutionary algorithms

The evolutionary algorithm used in this project (BRKGA) can be found inside `evolutionary_algorithms`, another git project.

## Scripts

Inside `scripts/` there are several useful scripts to perform different types of experiments.

## Test data

Test data can be found inside the `proteins/` directory. It contains data for several proteins. The data present for most of them is: fasta encoding (`fasta`), Rosetta fragments (`frag3` and `frag9`), PDB structure (`pdb`), secondary structure prediction (`ss_*`), and contact map prediction (`con_*`). Some proteins may contain extra info on `info` (i.e. the valid residue chain range that can be considered).

## Test results

Some test results can be found in `statistics`.