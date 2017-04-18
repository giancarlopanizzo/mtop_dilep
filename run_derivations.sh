#! /bin/sh
# run with parameter the path of this project ...
#  gSystem->Load("libExRootAnalysis.so");

filetoderive=$2
root -l -b <<- EOF
  gSystem->Load("libDelphes.so");
  .X ${1}/mtop_dilep/Derivations.C("""$filetoderive""");
  .q
EOF

