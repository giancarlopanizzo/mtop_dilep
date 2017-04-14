#! /bin/sh
# run with parameter the path of this project ...
#  gSystem->Load("libExRootAnalysis.so");

root -l -b <<- EOF
  gSystem->Load("libDelphes.so");
  .X ${1}/mtop_dilep/DelphesMoments.C();
  .q
EOF

