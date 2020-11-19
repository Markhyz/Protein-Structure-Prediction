#!/bin/bash

mv ./build/$1 $ROSETTA_BUILD_DIR/$1 2>/dev/null
./$ROSETTA_BUILD_DIR/$1 $2 $3 $4 $5