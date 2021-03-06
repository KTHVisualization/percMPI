#!/bin/bash -l
# The -l above is required to get the full environment with modules

build_type=$1

if [[ $build_type == "parallel" || $build_type == "multiple" ]]
then
    mkdir -p build
    cd build
    cmake -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS='-O3 -DCOMMUNICATION' ..
    make -j 8
elif [[ $build_type == "parallel-collectives" || $build_type == "multiple-collectives" ]]
then
    mkdir -p build_c
    cd build_c
    cmake -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS='-O3 -DCOMMUNICATION -DCOLLECTIVES' ..
    make -j 8
elif [[ $build_type == "sequential" || $build_type == "single" ]]
then
    mkdir -p build_s
    cd build_s
    cmake -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS='-O3 -DSINGLENODE' ..
    make -j 8
elif [[ $build_type == "sequentialcomm" ]]
then
    mkdir -p build_sc
    cd build_sc
    cmake -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS='-O3 -DCOMMUNICATION -DSINGLENODE' ..
    make -j 8

fi
