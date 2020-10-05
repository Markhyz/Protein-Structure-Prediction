#!/bin/bash

mv ./build/main $ROSETTA_BUILD_DIR/main 2>/dev/null
mpiexec --host localhost:$1 ./$ROSETTA_BUILD_DIR/main $2 $3 $4