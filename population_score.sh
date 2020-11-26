#!/bin/bash

mv ./build/population_score $ROSETTA_BUILD_DIR/population_score 2>/dev/null
./$ROSETTA_BUILD_DIR/population_score $1 $2 $3 $4 $5