#!/bin/bash

test_type=$1
account=$2
host=`uname -a`
user_email=`whoami`@kth.se

if [[ $host == *beskow* ]]
then
    num_procs_node=32
elif [[ $host == *tegner* ]]
then
    num_procs_node=24
else
    num_procs_node=4
fi

function makeTest
{
# test identifier
    test_name=$1
    datapath=$2
    rmsFile=$3 
    dataSizeX=$4 
    dataSizeY=$5
    dataSizeZ=$6
    timeStep=$7
    totalSizeX=$8
    totalSizeY=$9
    totalSizeZ=${10}
    blockSizeX=${11}
    blockSizeY=${12} 
    blockSizeZ=${13}
    hMin=${14} 
    hMax=${15} 
    hSamples=${16} 
    num_procs=${17}
# One global node + given number, make sure to round up
    num_nodes=$(((($num_procs + 1) + ($num_procs_node-1)) / $num_procs_node))
    job_name=$test_name"_n"$num_procs
    output_file="./"$job_name".sh"
    
    echo -e "#!/bin/bash -l\n# The -l above is required to get the full environment with modules\n" > $output_file
    echo -e "# Set the name of the script\n#SBATCH -J $job_name\n" >> $output_file
    echo -e "# Set the allocation to be charged for this job\n#SBATCH -A $account\n" >> $output_file
    echo -e "# Set the allocated time for the job\n#SBATCH -t 00:10:00\n" >> $output_file
    echo -e "# Set the node type to Haswell nodes only\n#SBATCH -C Haswell\n" >> $output_file
    echo -e "# Set the number of nodes\n#SBATCH --nodes=$num_nodes\n" >> $output_file
    echo -e "# Set the number of MPI processes\n#SBATCH -n $num_procs\n" >> $output_file
    echo -e "# Set the number of MPI processes per node\n#SBATCH --ntasks-per-node=$num_procs_node\n" >> $output_file
    echo -e "# Set the e-mail preferences for the user\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=$user_email\n" >> $output_file
    
    echo -e "echo -e \"Running job \\\\\"$job_name\\\\\"\\\\n---------------------\"" >> $output_file
    
    # Run one initial test (i.e., cold) and then 1 tests
    for test in `seq 1 1 2`
    do
        echo -e "echo -n \"    Performing test #$test... \"" >> $output_file
        echo -e "aprun -n $(($num_procs+1)) ./PercMPI --dataPath $datapath --rmsFile $rmsFile"\
            "--dataSize $dataSizeX $dataSizeY $dataSizeZ"\
	    "--timeStep $timeStep --totalSize $totalSizeX $totalSizeY $totalSizeZ"\
            "--blockSize $blockSizeX $blockSizeY $blockSizeZ"\
            "--hMin $hMin --hMax $hMax --hSamples $hSamples"\
            "--computeMode 0 --outputMode 1"\ >> $output_file
        echo -e "echo \"Done.\"" >> $output_file
    done
    
    echo -e "echo -e \"---------------------\\\\n\"" >> $output_file
    chmod a+x $output_file
}

# Make the strong-scaling tests for the 180-duct
if [[ $test_type == 1 ]]
then
    test_name="strong_duct180"
    datapath="../../Data/P3"
    rmsFile="uv_000"
    dataSizeX=193 
    dataSizeY=194
    dataSizeZ=1000
    timeStep=1
    totalSizeX=193
    totalSizeY=194
    totalSizeZ=1000
    hMin=0.0
    hMax=2.0 
    hSamples=1001
    for num_procs in 40 #512 256 128 32 16 8
    do
        blockSizeX=97
        blockSizeY=97 
        blockSizeZ=100
        makeTest $test_name $datapath $rmsFile $dataSizeX $dataSizeY $dataSizeZ $timeStep $totalSizeX $totalSizeY $totalSizeZ $blockSizeX $blockSizeY $blockSizeZ $hMin $hMax $hSamples $num_procs
    done
# Make the weak-scaling tests for the 180-duct
elif [[ $test_type == 2 ]]
then
    test_name="baseline_duct180"

fi

