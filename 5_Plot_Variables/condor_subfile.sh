Universe       = vanilla
# files will be copied back to this dir
Initialdir     = .
#tell condor where my grid certificate it
#x509userproxy=/tmp/x509up_u2142
# run my script
Executable     = condor_run.sh
+AccountingGroup = "group_cmshi.ilaflott"
#+IsMadgraph = 1
#Arguments      = $startfile $endfile $destination \$(Process)
# input files. in this case, there are none.
Input          = /dev/null
# log files
Error          = condor.stderr
Output         = condor.output
Log            = condor.log
# get the environment (path, etc.)
Getenv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
# should write all output & logs to a local directory
# and then transfer it back to Initialdir on completion
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
#specify any extra input files (for example, an orcarc file)
#transfer_input_files = $filelist
Queue
