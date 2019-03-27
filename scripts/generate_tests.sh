#!/bin/bash

test_type=$1
account=$2
host=`uname -a`
user_email=`whoami`@kth.se

if [[ $# < 2 ]]
then 
    echo "Not enough arguments. (Need Test Type and Account)"
fi


if [[ $host == *beskow* ]]
then
    num_procs_node=32
elif [[ $host == *tegner* ]]
then
    num_procs_node=24
else
    num_procs_node=32
fi

function makeTestDuct
{
# test identifier
    test_name=$1
    datapath="../Data/P3"
    rmsFile="uv_000"
    dataSizeX=193 
    dataSizeY=194
    dataSizeZ=1000
    timeStep=1
    totalSizeX=$2
    totalSizeY=$3
    totalSizeZ=$4
    blockSizeX=$5
    blockSizeY=$6 
    blockSizeZ=$7
    hMin=$8 
    hMax=$9 
    hSamples=${10} 
    num_procs=${11}
    build_dir=${12}
    # One global node + given number, make sure to round up
    if [[ $build_dir == "build_s" ]] 
    then 
        num_nodes=1
        # Incremented later by 1
        num_procs=0
    else
    num_nodes=$(((($num_procs + 1) + ($num_procs_node-1)) / $num_procs_node))
    fi
    job_name=$test_name"_n"$num_procs"_h"$hSamples
    num_procs_used=$(($num_procs+1))
    # evenly distribute the nodes
    if [[ $num_nodes == 1 ]]
    then
        num_procs_node_used=$num_procs_used
    else 
        num_procs_node_used=$((($num_procs_used + ($num_nodes-1)) / $num_nodes))
    fi
    output_file="./"$test_name"/"$job_name".sh"
    mkdir -p $test_name
    echo -e "#!/bin/bash -l\n# The -l above is required to get the full environment with modules\n" > $output_file
    echo -e "# Set the name of the script\n#SBATCH -J $job_name\n" >> $output_file
    echo -e "# Set the allocation to be charged for this job\n#SBATCH -A $account\n" >> $output_file
    echo -e "# Set the allocated time for the job\n#SBATCH -t 01:00:00\n" >> $output_file
    echo -e "# Set the node type to Haswell nodes only\n#SBATCH -C Haswell\n" >> $output_file
    echo -e "# Set the number of nodes\n#SBATCH --nodes=$num_nodes\n" >> $output_file
    echo -e "# Set the number of MPI processes\n#SBATCH -n $num_procs_used\n" >> $output_file
    echo -e "# Set the number of MPI processes per node\n#SBATCH --ntasks-per-node=$num_procs_node_used\n" >> $output_file
    echo -e "# Set the e-mail preferences for the user\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=$user_email\n" >> $output_file
    
    echo -e "echo -e \"Running job \\\\\"$job_name\\\\\"\\\\n---------------------\"" >> $output_file
    
    echo -e "echo -n \"    Performing initial test... \"" >> $output_file
        echo -e "aprun -n $num_procs_used -N $num_procs_node_used ./$build_dir/PercMPI --dataPath $datapath --rmsFile $rmsFile"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	        "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--inputMode 0 --computeMode 0 --outputMode 2 --outputPrefix $test_name"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    # Run 10 tests taking only timings
    for test in `seq 1 1 10`
    do
        echo -e "echo -n \"    Performing test #$test... \"" >> $output_file
        echo -e "aprun -n $num_procs_used -N $num_procs_node_used ./$build_dir/PercMPI --dataPath $datapath --rmsFile $rmsFile"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	        "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--inputMode 0 --computeMode 0 --outputMode 1 --outputPrefix $test_name"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    done
    
    echo -e "echo -e \"---------------------\\\\n\"" >> $output_file
    chmod a+x $output_file
}

function makeTestIso512
{
    # test identifier
    test_name=$1
    datapath=../Data/Iso4096/iso512x512x512.raw
    avgValue=0.0
    rmsValue=81.652634 
    dataSizeX=512 
    dataSizeY=512
    dataSizeZ=512
    timeStep=1
    totalSizeX=$2
    totalSizeY=$3
    totalSizeZ=$4
    blockSizeX=$5
    blockSizeY=$6 
    blockSizeZ=$7
    hMin=$8 
    hMax=$9 
    hSamples=${10} 
    num_procs=${11}
    build_dir=${12}
    # One global node + given number, make sure to round up
    if [[ $build_dir == "build_s" ]] 
    then 
        num_nodes=1
        # Incremented later by 1
        num_procs=0
    else
    num_nodes=$(((($num_procs + 1) + ($num_procs_node-1)) / $num_procs_node))
    fi
    job_name=$test_name"_n"$num_procs"_h"$hSamples
    num_procs_used=$(($num_procs+1))
    # evenly distribute the nodes
    if [[ $num_nodes == 1 ]]
    then
        num_procs_node_used=$num_procs_used
    else 
        num_procs_node_used=$((($num_procs_used + ($num_nodes-1)) / $num_nodes))
    fi
    output_file="./"$test_name"/"$job_name".sh"
    mkdir -p $test_name
    echo -e "#!/bin/bash -l\n# The -l above is required to get the full environment with modules\n" > $output_file
    echo -e "# Set the name of the script\n#SBATCH -J $job_name\n" >> $output_file
    echo -e "# Set the allocation to be charged for this job\n#SBATCH -A $account\n" >> $output_file
    echo -e "# Set the allocated time for the job\n#SBATCH -t 02:00:00\n" >> $output_file
    echo -e "# Set the node type to Haswell nodes only\n#SBATCH -C Haswell\n" >> $output_file
    echo -e "# Set the number of nodes\n#SBATCH --nodes=$num_nodes\n" >> $output_file
    echo -e "# Set the number of MPI processes\n#SBATCH -n $num_procs_used\n" >> $output_file
    echo -e "# Set the number of MPI processes per node\n#SBATCH --ntasks-per-node=$num_procs_node_used\n" >> $output_file
    echo -e "# Set the e-mail preferences for the user\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=$user_email\n" >> $output_file
    
    echo -e "echo -e \"Running job \\\\\"$job_name\\\\\"\\\\n---------------------\"" >> $output_file
    
    echo -e "echo -n \"    Performing initial test... \"" >> $output_file
        echo -e "aprun -n $num_procs_used -N $num_procs_node_used ./$build_dir/PercMPI --dataPath $datapath --avgValue $avgValue --rmsValue $rmsValue"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	        "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--inputMode 3 --computeMode 0 --outputMode 2 --outputPrefix $test_name"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    # Run 10 tests taking only timings
    for test in `seq 1 1 10`
    do
        echo -e "echo -n \"    Performing test #$test... \"" >> $output_file
        echo -e "aprun -n $num_procs_used -N $num_procs_node_used ./$build_dir/PercMPI --dataPath $datapath --avgValue $avgValue --rmsValue $rmsValue"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	        "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--inputMode 3 --computeMode 0 --outputMode 1 --outputPrefix $test_name"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    done
    
    echo -e "echo -e \"---------------------\\\\n\"" >> $output_file
    chmod a+x $output_file
}

function makeTestIso4096
{
    # test identifier
    test_name=$1
    datapath=../Data/Iso4096/iso4096.raw
    avgValue=0.0
    rmsValue=81.652634 
    dataSizeX=4089
    dataSizeY=4089
    dataSizeZ=4089
    timeStep=1
    totalSizeX=$2
    totalSizeY=$3
    totalSizeZ=$4
    blockSizeX=$5
    blockSizeY=$6 
    blockSizeZ=$7
    hMin=$8 
    hMax=$9 
    hSamples=${10} 
    num_procs=${11}
    build_dir=${12}
    # One global node + given number, make sure to round up
    if [[ $build_dir == "build_s" ]] 
    then 
        num_nodes=1
        # Incremented later by 1
        num_procs=0
    else
    num_nodes=$(((($num_procs + 1) + ($num_procs_node-1)) / $num_procs_node))
    fi
    job_name=$test_name"_n"$num_procs"_h"$hSamples
    num_procs_used=$(($num_procs+1))
    # evenly distribute the nodes
    if [[ $num_nodes == 1 ]]
    then
        num_procs_node_used=$num_procs_used
    else 
        num_procs_node_used=$((($num_procs_used + ($num_nodes-1)) / $num_nodes))
    fi
    output_file="./"$test_name"/"$job_name".sh"
    mkdir -p $test_name
    echo -e "#!/bin/bash -l\n# The -l above is required to get the full environment with modules\n" > $output_file
    echo -e "# Set the name of the script\n#SBATCH -J $job_name\n" >> $output_file
    echo -e "# Set the allocation to be charged for this job\n#SBATCH -A $account\n" >> $output_file
    echo -e "# Set the allocated time for the job\n#SBATCH -t 02:00:00\n" >> $output_file
    echo -e "# Set the node type to Haswell nodes only\n#SBATCH -C Haswell\n" >> $output_file
    echo -e "# Set the number of nodes\n#SBATCH --nodes=$num_nodes\n" >> $output_file
    echo -e "# Set the number of MPI processes\n#SBATCH -n $num_procs_used\n" >> $output_file
    echo -e "# Set the number of MPI processes per node\n#SBATCH --ntasks-per-node=$num_procs_node_used\n" >> $output_file
    echo -e "# Set the e-mail preferences for the user\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=$user_email\n" >> $output_file
    
    echo -e "echo -e \"Running job \\\\\"$job_name\\\\\"\\\\n---------------------\"" >> $output_file
    
    echo -e "echo -n \"    Performing initial test... \"" >> $output_file
        echo -e "aprun -n $num_procs_used -N $num_procs_node_used ./$build_dir/PercMPI --dataPath $datapath --avgValue $avgValue --rmsValue $rmsValue"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	        "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--inputMode 3 --computeMode 0 --outputMode 2 --outputPrefix $test_name"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    # Run 10 tests taking only timings
    for test in `seq 1 1 0` # Todo change this for more samples
    do
        echo -e "echo -n \"    Performing test #$test... \"" >> $output_file
        echo -e "aprun -n $num_procs_used -N $num_procs_node_used ./$build_dir/PercMPI --dataPath $datapath --avgValue $avgValue --rmsValue $rmsValue"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	        "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--inputMode 3 --computeMode 0 --outputMode 1 --outputPrefix $test_name"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    done
    
    echo -e "echo -e \"---------------------\\\\n\"" >> $output_file
    chmod a+x $output_file
}

function makeTestMergeIso512
{
    # test identifier
    test_name=$1
    datapath=../Data/Iso4096/iso512x512x512.vti
    executablepath=./../FTC/build/standalone/FTCTree/cmd/ftcTreeCmd
    num_procs=$2
    job_name=$test_name"_n"$num_procs
    output_file="./"$test_name"/"$job_name".sh"
    mkdir -p $test_name
    echo -e "#!/bin/bash -l\n# The -l above is required to get the full environment with modules\n" > $output_file
    echo -e "# Set the name of the script\n#SBATCH -J $job_name\n" >> $output_file
    echo -e "# Set the allocation to be charged for this job\n#SBATCH -A $account\n" >> $output_file
    echo -e "# Set the allocated time for the job\n#SBATCH -t 00:20:00\n" >> $output_file
    echo -e "# Set the node type to Haswell nodes only\n#SBATCH -C Haswell\n" >> $output_file
    echo -e "# Set the number of nodes\n#SBATCH --nodes=1\n" >> $output_file
    echo -e "# Set the number of MPI processes\n#SBATCH -n 1\n" >> $output_file
    echo -e "# Set the number of MPI processes per node\n#SBATCH --ntasks-per-node=1\n" >> $output_file
    echo -e "# Set the e-mail preferences for the user\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=$user_email\n" >> $output_file
    echo -e "# Set the number of OpenMP processes\n#SBATCH -c $num_procs\n"  >> $output_file
    echo -e "# Set outfile\n#SBATCH -o $job_name.out$\n"  >> $output_file
    echo -e "echo -e \"Running job \\\\\"$job_name\\\\\"\\\\n---------------------\"" >> $output_file
    
    echo -e "export OMP_NUM_THREADS=$num_procs"  >> $output_file

    echo -e "echo -n \"    Performing initial test... \"" >> $output_file
    echo -e "aprun -n 1 -N 1 -d $num_procs $executablepath -g $datapath -d 3 -t $num_procs -f 0 -T 0" >> $output_file
    echo -e "echo \"Done.\"" >> $output_file
    # Run 10 tests taking only timings
    for test in `seq 1 1 9` # Change for 
    do
        echo -e "echo -n \"    Performing test #$test... \"" >> $output_file
        echo -e "aprun -n 1 -N 1 -d $num_procs $executablepath -g $datapath -d 3 -t $num_procs -f 0 -T 0" >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    done
    
    echo -e "echo -e \"---------------------\\\\n\"" >> $output_file
    chmod a+x $output_file
}


# Make the strong-scaling tests for the 180-duct
if [[ $test_type == 1 ]]
then
    test_name="strong_duct180"
    echo $test_name
    totalSizeX=193
    totalSizeY=194
    totalSizeZ=1000
    hMin=0.0
    hMax=2.0 
    hSamples=1000
    blockSizeX=193
    blockSizeY=194 
    blockSizeZ=1000
    num_procs=1
    for tests in `seq 1 1 9`
    do
        num_procs=$(($num_procs*2))
        if [[ $blockSizeX -ge $blockSizeY && $blockSizeX -ge $blockSizeZ ]]
        then
            blockSizeX=$((($blockSizeX+1)/2))
        elif [[ $blockSizeY -ge $blockSizeX && $blockSizeY -ge $blockSizeZ ]]
        then
            blockSizeY=$((($blockSizeY+1)/2))
        else
            blockSizeZ=$((($blockSizeZ+1)/2))
        fi
        echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
        for hSamples in 100 1000
        do
        makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
        done
    done
# Make the weak-scaling tests for the 180-duct
elif [[ $test_type == 2 ]]
then
    test_name="baseline_duct180"
    echo $test_name
    totalSizeX=193
    totalSizeY=194
    totalSizeZ=1000
    hMin=0.0
    hMax=2.0 
    hSamples=1
    blockSizeX=193
    blockSizeY=194 
    blockSizeZ=1000
    num_procs=1
    for tests in `seq 1 1 5`
    do
        hSamples=$(($hSamples*10))
        echo "Generating baseline test with "$hSamples" samples."
        makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build_s"
    done
# Make the weak-scaling tests for the 180-duct
elif [[ $test_type == 3 ]]
then
    test_name="weak_duct180"
    echo $test_name
    hMin=0.0
    hMax=2.0 
    blockSizeX=48
    blockSizeY=48 
    blockSizeZ=31
    totalSizeX=$blockSizeX
    totalSizeY=$blockSizeY
    totalSizeZ=$blockSizeZ
    num_procs=1
    # Test for the single block
    for hSamples in 100 1000
    do
    echo "Generating weak scaling test for "$num_procs" processes with total ("$totalSizeX $totalSizeY $totalSizeZ")."
    makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build_s"
    done
    num_procs=1
    for tests in `seq 1 1 9`
    do
        num_procs=$(($num_procs*2))
        if [[ $(($totalSizeX*2)) -le $dataSizeX && $totalSizeX -le $totalSizeY && $totalSizeX -le $totalSizeZ ]]
        then
            totalSizeX=$(($totalSizeX*2))
        elif [[ $(($totalSizeY*2)) -le $dataSizeY && $totalSizeY -le $totalSizeX && $totalSizeY -le $totalSizeZ ]]
        then
            totalSizeY=$(($totalSizeY*2))
        else
            totalSizeZ=$(($totalSizeZ*2))
        fi
        echo "Generating weak scaling test for "$num_procs" processes with total ("$totalSizeX $totalSizeY $totalSizeZ")."
        for hSamples in 100 1000
        do
        makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
        done
    done
# Make the weak-scaling tests for the 180-duct
elif [[ $test_type == 4 ]]
then
    test_name="weakh_duct180"
    echo $test_name
    hMin=0.0
    hMax=2.0 
    blockSizeX=193
    blockSizeY=97 
    blockSizeZ=125
    totalSizeX=193
    totalSizeY=194
    totalSizeZ=1000
    num_procs=16
    hSamples=1
    for tests in `seq 1 1 5`
    do
        hSamples=$((10*$hSamples))
        echo "Generating weak scaling test with "$hSamples" samples."
        makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    done
# Make the strong-scaling tests for the iso512
elif [[ $test_type == 5 ]]
then
    test_name="strong_iso512"
    echo $test_name
    hMin=0.0
    hMax=4.0 
    hSamples=1000
    blockSizeX=512
    blockSizeY=512 
    blockSizeZ=512
    totalSizeX=512
    totalSizeY=512
    totalSizeZ=512
    num_procs=1
    echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
    makeTestIso512 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build_s"
    num_procs=1
    # Create tests for powers of two
    for tests in `seq 1 1 9`
    do
        num_procs=$(($num_procs*2))
        if [[ $blockSizeX -ge $blockSizeY && $blockSizeX -ge $blockSizeZ ]]
        then
            blockSizeX=$((($blockSizeX+1)/2))
        elif [[ $blockSizeY -ge $blockSizeX && $blockSizeY -ge $blockSizeZ ]]
        then
            blockSizeY=$((($blockSizeY+1)/2))
        else
            blockSizeZ=$((($blockSizeZ+1)/2))
        fi
        echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
        makeTestIso512 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    done
    # Create tests for not powers of two (to compare against merge tree
    num_procs=6 #2*3
    blockSizeX=256
    blockSizeY=171 
    blockSizeZ=512
    echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
        makeTestIso512 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    blockSizeX=256
    blockSizeY=103 
    blockSizeZ=512
    num_procs=10 #2*5
    echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
        makeTestIso512 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"    
    blockSizeX=256
    blockSizeY=256 
    blockSizeZ=171
    num_procs=12 #2*2*3
    echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
        makeTestIso512 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    num_procs=14 #2*7
    blockSizeX=256
    blockSizeY=74 
    blockSizeZ=512
    echo "Generating strong scaling test for "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
        makeTestIso512 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
# Make the weak-scaling tests for the 180-duct
elif [[ $test_type == 6 ]]
then
    test_name="weak_iso512"
    echo $test_name
    hMin=0.0
    hMax=4.0 
    blockSizeX=128
    blockSizeY=128 
    blockSizeZ=128
    totalSizeX=$blockSizeX
    totalSizeY=$blockSizeY
    totalSizeZ=$blockSizeZ
    num_procs=1
    hsamples=1000
    # Test for the single block
    echo "Generating weak scaling test for "$num_procs" processes with total ("$totalSizeX $totalSizeY $totalSizeZ")."
    makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build_s"
    num_procs=1
    for tests in `seq 1 1 9`
    do
        num_procs=$(($num_procs*2))
        if [[ $(($totalSizeX*2)) -le 512 && $totalSizeX -le $totalSizeY && $totalSizeX -le $totalSizeZ ]]
        then
            totalSizeX=$(($totalSizeX*2))
        elif [[ $(($totalSizeY*2)) -le 512 && $totalSizeY -le $totalSizeX && $totalSizeY -le $totalSizeZ ]]
        then
            totalSizeY=$(($totalSizeY*2))
        else
            totalSizeZ=$(($totalSizeZ*2))
        fi
        echo "Generating weak scaling test for "$num_procs" processes with total ("$totalSizeX $totalSizeY $totalSizeZ")."
        makeTestDuct $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    done
# Make the strong-scaling tests for the iso512
elif [[ $test_type == 7 ]]
then
    test_name="strong_iso4096"
    echo $test_name
    hMin=0.0
    hMax=4.0 
    hSamples=1000
    totalSizeX=4089
    totalSizeY=4089
    totalSizeZ=4089
    nodesX=4
    nodesY=4
    nodesZ=2
    num_procs=$(($nodesX*$nodesY*$nodesZ*4*4*2))
    blockSizeX=$((($totalSizeX+($nodesX*4-1))/($nodesX*4)))
    blockSizeY=$((($totalSizeY+($nodesY*4-1))/($nodesY*4)))
    blockSizeZ=$((($totalSizeZ+($nodesZ*2-1))/($nodesZ*2)))
    echo "Generating strong scaling test for "$(($nodesX*$nodesY*$nodesZ))" nodes, "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
    makeTestIso4096 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    nodesX=4
    nodesY=3
    nodesZ=3
    num_procs=$(($nodesX*$nodesY*$nodesZ*4*4*2))
    blockSizeX=$((($totalSizeX+($nodesX*4-1))/($nodesX*4)))
    blockSizeY=$((($totalSizeY+($nodesY*4-1))/($nodesY*4)))
    blockSizeZ=$((($totalSizeZ+($nodesZ*2-1))/($nodesZ*2)))
    echo "Generating strong scaling test for "$(($nodesX*$nodesY*$nodesZ))" nodes, "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
    makeTestIso4096 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    nodesX=4
    nodesY=4
    nodesZ=3
    num_procs=$(($nodesX*$nodesY*$nodesZ*4*4*2))
    blockSizeX=$((($totalSizeX+($nodesX*4-1))/($nodesX*4)))
    blockSizeY=$((($totalSizeY+($nodesY*4-1))/($nodesY*4)))
    blockSizeZ=$((($totalSizeZ+($nodesZ*2-1))/($nodesZ*2)))
    echo "Generating strong scaling test for "$(($nodesX*$nodesY*$nodesZ))" nodes, "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
    makeTestIso4096 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
    nodesX=4
    nodesY=4
    nodesZ=4
    num_procs=$(($nodesX*$nodesY*$nodesZ*4*4*2))
    blockSizeX=$((($totalSizeX+($nodesX*4-1))/($nodesX*4)))
    blockSizeY=$((($totalSizeY+($nodesY*4-1))/($nodesY*4)))
    blockSizeZ=$((($totalSizeZ+($nodesZ*2-1))/($nodesZ*2)))
    echo "Generating strong scaling test for "$(($nodesX*$nodesY*$nodesZ))" nodes, "$num_procs" processes with blocksize ("$blockSizeX $blockSizeY $blockSizeZ")."
    makeTestIso4096 $test_name $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs "build"
elif [[ $test_type == 8 ]]
then
    test_name="merge_iso512"
    echo $test_name
    echo "Generating scaling test for 1 for merge tree computation."
    makeTestMergeIso512 $test_name 1
    for num_procs in `seq 2 2 16`
    do
    echo "Generating scaling test for "$num_procs" for merge tree computation."
    makeTestMergeIso512 $test_name $num_procs
    done
fi

