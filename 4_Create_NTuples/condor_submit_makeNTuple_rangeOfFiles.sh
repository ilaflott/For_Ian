#!/bin/bash

# run this with $ source ppNTuples_condor_submit.sh 
#if arguments use $ source ppNTuples_condor_submit.sh $arg1 $arg2 ... etc.
#echo "inputs are \$flavor \$BeginJob \$EndJob \$JobSplittng"
#echo "to run one job let BeginJob=EndJob"

job=0
flavor=$1
BeginJob=$2
EndJob=$3
JobSplitting=$4

#echo "$NJobSplitting being submitted"
#echo "makeNTuple Jobs Being Submitted"

HadoopDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples"
NTupleDir=""
NTupleFileName=""
FinalNTupleFileName=""

if [ $flavor -eq 0 ]; then
#    echo "Data Job Being Submitted"
    NTupleDir="${HadoopDir}/Data_redo_11.1"
    NTupleFileName="${NTupleDir}/data_NTuple_11.1_reForest"
#    NTupleDir="{$HadoopDir}/Data_MaybeBadAndOldFiles_11.1"
#    NTupleFileName="${NTupleDir}/data_NTuple_11.1__"
fi
if [ $flavor -eq 1 ]; then
#    echo "QCDJet Job Submitted"
    NTupleDir="${HadoopDir}/QCDJets_Official_noVsJets"
    NTupleFileName="${NTupleDir}/QCDJets_NTuple_noVsJets_11.1_"
fi
if [ $flavor -eq 2 ]; then
#    echo "BJet Job Being Submitted"
    NTupleDir="${HadoopDir}/BJets_OfficialLowPtAndPrivateGenHighPT_11.1"
#    NTupleDir="${HadoopDir}/BJets_addStat_11.1"
    NTupleFileName="${NTupleDir}/BJets_NTuple_addStat_11.1_"
fi
if [ $flavor -eq 3 ]; then
#    echo "CJet Job Being Submitted"
    NTupleDir="${HadoopDir}/CJets_OfficialLowPtAndPrivateGenHighPT_11.1"
#    NTupleDir="${HadoopDir}/CJets_addStat_11.1"
    NTupleFileName="${NTupleDir}/CJets_NTuple_addStat_11.1_"
fi

JobNum=$BeginJob

while [ $JobNum -le $EndJob ]
do
#  echo "submitting Job number $JobNum"
    
    cat > subfile <<EOF
Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
# files will be copied back to this dir
#Initialdir     = .
# tell condor where my grid certificate is if needed
#x509userproxy=/tmp/x509up_u2142
# run my script
Executable     = condor_run_makeNTuple.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $flavor $JobNum $JobSplitting
#input files. in this case, there are none.
Input          = /dev/null
# log files
Error          = Condor_logs/makeNTuple_Flav_${flavor}_p${JobNum}_of_${JobSplitting}.err
Output         = Condor_logs/makeNTuple_Flav_${flavor}_p${JobNum}_of_${JobSplitting}.out
Log            = Condor_logs/makeNTuple_Flav_${flavor}_p${JobNum}_of_${JobSplitting}.log
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
    FinalNTupleFileName="${NTupleFileName}_${JobNum}_of_${JobSplitting}.root"
    actualsize=$(wc -c <"${FinalNTupleFileName}")

#    echo "actualsize is ${actualsize}"#debug
    if [ ! -e $FinalNTupleFileName ]; then
	echo "file ${FinalNTupleFileName} does not exist!"
	echo "submitting subfile..."
	condor_submit  subfile
	sleep 60s
    elif [ $actualsize -le 8500 ]; then
	echo "file ${FinalNTupleFileName} does exist but isn't big enough!"
	echo "submitting subfile..."
	condor_submit  subfile
	sleep 60s
    else
#	echo "file ${FinalNTupleFileName} exists and is big enough! DON'T SUBMIT"
	JobNum=$(($JobNum + 1))
	continue
    fi
    JobNum=$(($JobNum + 1))
done