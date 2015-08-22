#!bin/bash
#this macro makes the tarball that the submit script points condor to

#compile the macro
root -l <<EOF
.L bTagNTuple.C++
.q
EOF

#remove tarball if already there so it can be remade
rm run_ppNTuples.tar

#tar [options] tarballname  input_files  macro-code+dependencies from compilation
tar -zcvf run_ppNTuples.tar filelists/*Jets_filelist*.txt filelists/ppMuon2013A_runForest_filelist.txt weight_info/*.txt bTagNTuple.C bTagNTuple_C.* 
