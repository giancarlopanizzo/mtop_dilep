#! /bin/sh
# run with parameter the path of this project ...
#  gSystem->Load("libExRootAnalysis.so");

combine=$2
fixmp=$3
lhaid=$4
root -l -b <<- EOF
  gSystem->Load("libDelphes.so");
  .X ${1}/mtop_dilep/DelphesMoments.C($combine,$lhaid,$fixmp);
  .q
EOF

