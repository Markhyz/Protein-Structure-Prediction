#!/bin/bash

mv ./build/population_score $ROSETTA_BUILD_DIR/population_score 2>/dev/null
./$ROSETTA_BUILD_DIR/population_score $@
