#! /bin/sh

root -l -b <<- EOF
  gSystem->Load("libExRootAnalysis.so");
  .X SimpleDraw.C();
  .q
EOF

