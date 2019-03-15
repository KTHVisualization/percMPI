#!/bin/bash

script=$1
path=$(dirname "$0")/jobid.tmp
prev_jobid=`cat $path 2> /dev/null`

if [[ $prev_jobid > 0 ]]
then
    echo -n "Queuing script \"$script\" after jobid=$prev_jobid... "
    job=`sbatch --dependency=afterany:$prev_jobid $script`
else
    echo -n "Queuing a new job with the script \"$script\"... "
    job=`sbatch $script`
fi

# Store the current jobid for future sbatch calls
jobid=`echo $job | awk '{print $4}'`
echo -n $jobid > $path

echo "Done."

