#! /bin/sh

root -l -b <<- EOF
  gSystem->Load("libExRootAnalysis.so");
  gSystem->Load("libDelphes.so");
  .X /home/gian/Public/CERN//Projects/mtop_dilep/DelphesMoments.C();
  .q
EOF

