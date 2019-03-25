#!/bin/bash

dir = $1

for job in $1*.sh 
do
    ./scripts/sbatch.sh $job
done 
