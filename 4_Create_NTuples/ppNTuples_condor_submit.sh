#!/bin/bash

# run this with $ source ppNTuples_condor_submit.sh 
#if arguments use $ source ppNTuples_condor_submit.sh $arg1 $arg2 ... etc.

#job=0
#flavor=4
#max=14
echo "inputs are \$job \$flavor \$maxFlavor"
job=$1
flavor=$2
max=$3
NJOBS=$(($(($max + 1)) - $flavor)) 

echo "$NJOBS being submitted"

    if [ $job -eq 0 ]; then
      echo "makeNTuple Jobs Being Submitted"
    fi
    if [ $job -eq 1 ]; then
      echo "MCCounts Jobs Being Submitted"
    fi

while [ $flavor -le $max ]
do
      if [ $flavor -eq 0 ]; then
        echo "Data Job Being Submitted"
      fi
      if [ $flavor -eq 1 ]; then
        echo "FULL QCDJet Job Being Submitted"
        echo "Are you sure you want to do this?!"
      fi
      if [ $flavor -eq 2 ]; then
        echo "BJet Being Submitted"
      fi
      if [ $flavor -eq 3 ]; then
        echo "CJet Being Submitted"
      fi
      if [ $flavor -ge 4 ]; then
	Piece=$(($flavor - 3))  
        echo "QCDJet $Piece Being Submitted"
      fi
  echo "submitting $flavor"
    
    cat > subfile <<EOF
Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
# files will be copied back to this dir
#Initialdir     = .
# tell condor where my grid certificate is if needed
#x509userproxy=/tmp/x509up_u2142
# run my script
Executable     = ppNTuples_condor_run.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $job $flavor
#input files. in this case, there are none.
Input          = /dev/null
# log files
Error          = Condor_logs/Condor_Job_${job}_Flav_${flavor}.err
Output         = Condor_logs/Condor_Job_${job}_Flav_${flavor}.out
Log            = Condor_logs/Condor_Job_${job}_Flav_${flavor}.log
# get the environment (path, variables set in this submit script etc.)
Getenv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
transfer_input_files = run_ppNTuples.tar
# should write all output & logs to a local directory
# and then transfer it back to Initialdir on completion
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
# specify any extra input files (for example, an orcarc file)
Queue
EOF
    #submit the job
    echo "submitting subfile..."
    condor_submit -name "bTNT_${job}_${flavor}" subfile
    echo "Condor Job name is bTNT_${job}_${flavor}"
    flavor=$(($flavor + 1))
done
