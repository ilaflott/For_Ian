#!/bin/bash

#untar the tarball fed to condor via the submit script
#echo "untarring tarball"
tar -zxvf run_ppNTuples.tar > /dev/null

#cmsenv
echo "cmsenv"
#export SCRAM_ARCH=slc6_amd64_gcc472
source /osg/app/cmssoft/cms/cmsset_default.sh

#something i'll leave here for now
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

#what compiler is being used
gcc --version

#debugging
echo "listing contents.."
ls

#arguments used by macro fed to runscript by submit script
echo "Processing..."
job=0
flavor=$1
JobNum=$2
NJobs=$3

###for redoing pieces of the analysis
#rm /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_weight_info/*.txt
#rm /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/*.root
echo "Running makeNTuple Job"

#additional instructions for specific job...
if [ $JobNum -eq 1 ]; then
    echo "running other instructions too.."
fi

NTupleDir=""
HadoopDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples"
if [ $flavor -eq 0 ]; then
#    NTupleDir="${HadoopDir}/Data_MaybeBadAndOldFiles_11.1"
    NTupleDir="${HadoopDir}/Data_redo_11.1"
#    NTupleDir="${HadoopDir}/TESTING"
fi
if [ $flavor -eq 1 ]; then
    NTupleDir="${HadoopDir}/QCDJets_Official_noVsJets"
fi
if [ $flavor -eq 2 ]; then
#    NTupleDir="${HadoopDir}/BJets_addStat_11.1"
#    NTupleDir="${HadoopDir}/BJets_OfficialLowPt_11.1
#    NTupleDir="${HadoopDir}/BJets_HighPtOnly_11.1
    NTupleDir="${HadoopDir}/BJets_OfficialLowPtAndPrivateGenHighPT_11.1"
fi
if [ $flavor -eq 3 ]; then
#    NTupleDir="${HadoopDir}/CJets_addStat_11.1"
#    NTupleDir="${HadoopDir}/CJets_OfficialLowPt_11.1
#    NTupleDir="${HadoopDir}/CJets_HighPtOnly_11.1
    NTupleDir="${HadoopDir}/CJets_OfficialLowPtAndPrivateGenHighPT_11.1"
fi

#run the code
root -b -l <<EOF
.x bTagNTuple.C+(${job}, ${flavor},${JobNum},${NJobs})
.q
EOF

echo "moving output to hadoop..."
mv *.root "${NTupleDir}"

echo "done!"