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
NTupleFileName=""
if [ $flavor -eq 0 ]; then
    echo "moving Data file to Data Folder..."
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/Data_MaybeBadAndOldFiles_11.1"
    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/Data_redo_11.1"
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/TESTING"
    NTupleFileName="${NTupleDir}/data_NTuple_11.1_reForest_${JobNum}_of_${NJobs}.root"
fi
if [ $flavor -eq 1 ]; then
    echo "moving QCDJets file to QCDJets Folder..."
    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/QCDJets_Official_noVsJets"
    NTupleFileNAme=""
fi
if [ $flavor -eq 2 ]; then
    echo "moving BJets file to BJets Folder..."
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/BJets_addStat_11.1"
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/BJets_OfficialLowPt_11.1
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/BJets_HighPtOnly_11.1
    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/BJets_OfficialLowPtAndPrivateGenHighPT_11.1"
fi
if [ $flavor -eq 3 ]; then
    echo "moving CJets file to CJets Folder..."
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/CJets_addStat_11.1"
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/CJets_OfficialLowPt_11.1
#    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/CJets_HighPtOnly_11.1
    NTupleDir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/CJets_OfficialLowPtAndPrivateGenHighPT_11.1"
fi

#if file exists...
if [ ! -e $NTupleFileName ]; then
#run the code
    root -b -l <<EOF
.x bTagNTuple.C+(${job}, ${flavor},${JobNum},${NJobs})
.q
EOF
mv *.root "${NTupleDir}"
else
#file exists, don't run the code!
    echo "file exists!"
fi

echo "done!"