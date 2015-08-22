#!/bin/bash

#untar the tarball fed to condor via the submit script
#echo "untarring tarball"
tar -zxvf run_ppNTuples.tar > /dev/null

#cmsenv
echo "cmsenv"
export SCRAM_ARCH=slc6_amd64_gcc472
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
job=$1
flavor=$2

echo "whoami"
whoami
echo "running hadoop permissions test..."
echo "attempting to make directory..."
mkdir /mnt/hadoop/cms/store/user/ilaflott/test
echo "attempting to touch"
touch /mnt/hadoop/cms/store/user/ilaflott/test.txt
echo "attempting a simple copy commmand"
cp filelists/QCDJets_noVsJets_filelist_1.txt /mnt/hadoop/cms/store/user/ilaflott

#execute the script like this...
#root -b -l <<EOF
#.x bTagNTuple.C+($job, $flavor)
#.q
#EOF
##...not like this
# root -b -q bTagNTuple_Original.C\+\($flavor\) 

echo "Done!"

