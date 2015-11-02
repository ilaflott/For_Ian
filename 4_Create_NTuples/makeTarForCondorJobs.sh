#!/bin/bash
#this macro makes the tarball that the submit script points condor to

#compile the macro
root -l <<EOF
.L bTagNTuple.C++
.q
EOF

#remove tarball if already there so it can be remade
rm run_ppNTuples.tar

#tar [options] tarballname  input_files  macro-code+dependencies from compilation
tar -zcvf run_ppNTuples.tar filelists/*.txt  bTagNTuple.C bTagNTuple_C.* 
