#!/bin/bash

# run this with $ source ppNTuples_condor_submit.sh 
#if arguments use $ source ppNTuples_condor_submit.sh $arg1 $arg2 ... etc.

echo "inputs are \$flavor \$BeginJob \$EndJob \$JobSplittng"
echo "to run one job let BeginJob=EndJob"
job=0
flavor=$1
BeginJob=$2
EndJob=$3
JobSplitting=$4
#jobSegment=$3

#echo "$NJobSplitting being submitted"
echo "makeNTuple Jobs Being Submitted"


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

JobNum=$BeginJob
echo '${JobNum}'
while [ $JobNum -le $EndJob ]
#while [ $JobNum -le 5 ]
do
    
    root -b -l <<EOF
.x bTagNTuple.C+(${job},${flavor},${JobNum},${JobSplitting})
.q
EOF
    
    JobNum=$(($JobNum + 1))
#    echo "sleeping for 40s..."
#    sleep 10s
    echo "moving file"
    mv *.root /net/hidsk0001/d00/scratch/ilaflott/pp_NTuples/oldStat_HF_Pieces
#    echo "sleeping for 40s..."
#    sleep 10s
done