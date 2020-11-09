#!/bin/bash

mv ./build/score_calc $ROSETTA_BUILD_DIR/score_calc 2>/dev/null
./$ROSETTA_BUILD_DIR/score_calc $1 $2