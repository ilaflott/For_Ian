#!/bin/sh

#export X509_USER_PROXY=/net/hisrv0001/home/rkunnawa/myproxy/
joblist=runTheseCondorJobs.list
nJobs=`wc -l < $joblist`
counter=1
echo "submitting $nJobs..."
while read LINE 
do
    # Condor run file
    echo $LINE
    cat > runfile <<EOF
export SCRAM_ARCH=slc5_amd64_gcc462
source /osg/app/cmssoft/cms/cmsset_default.sh
cd /net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src
eval `scramv1 runtime -sh`
cd /net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/5_Plot_Variables/
echo "Processing..."
$LINE
echo "Done!"
EOF
    cat > subfile <<EOF
Universe       = vanilla
Initialdir     = .
#tell condor where my grid certificate is if it needs it
#x509userproxy=/tmp/x509up_u2142
Executable     = runfile
+AccountingGroup = "group_cmshi.ilaflott"
Input          = /dev/null
Error          = condor.stderr
Output         = condor.output
Log            = condor.log
Getenv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
# should write all output & logs to a local directory
# and then transfer it back to Initialdir on completion
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Queue
EOF
    # submit the job
    echo "submitting job# $counter... " 
    condor_submit subfile
    counter=$(($counter+1))
done<$joblist