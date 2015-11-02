#!/bin/bash

# run this with $ source ppNTuples_condor_submit.sh 
#if arguments use $ source ppNTuples_condor_submit.sh $arg1 $arg2 ... etc.

#job=0
#flavor=4
#max=14
echo "inputs are \$flavor"
flavor=$1
#jobSegment=$3
#NJobs=$3

echo "hadd being submitted"

if [ $flavor -eq 0 ]; then
    echo "Data Job Being Submitted"
fi
if [ $flavor -eq 1 ]; then
    echo "QCDJet Job Submitted"
fi
if [ $flavor -eq 2 ]; then
    echo "BJet Job Being Submitted"
fi
if [ $flavor -eq 3 ]; then
    echo "CJet Job Being Submitted"
fi

cat > subfile <<EOF
Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
# files will be copied back to this dir
#Initialdir     = .
# tell condor where my grid certificate is if needed
#x509userproxy=/tmp/x509up_u2142
# run my script
Executable     = condor_run_haddJob.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $flavor
#input files. in this case, there are none.
Input          = /dev/null
# log files
Error          = Condor_logs/haddJob_Flav_${flavor}.err
Output         = Condor_logs/haddJob_Flav_${flavor}.out
Log            = Condor_logs/haddJob_Flav_${flavor}.log
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
    #bTNT->bTagNTuple
echo "Condor Job name is bTNT_haddJob_${flavor}"
condor_submit -name "bTNT_haddJob_${flavor}" subfile

done