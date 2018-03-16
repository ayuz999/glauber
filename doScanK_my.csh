#!/bin/csh
set nppbin  = "1"
set nppmin  = $1
set nppmax  = $1

set xbin    = "1"
set xmin    = $2
set xmax    = $2

set kbin    = "1"
set kmin    = $3
set kmax    = $3

set efficiency = $4 # for multiplicity dependent efficiency

starver SL16a

set nevents  = "1000000" # 1M
set real     = "refmult3.root"
#set real     = "npart.root"
set mc       = "ncoll_npart.root"
set multCut  = "100" #for low multiplicity cut
#set efficiency = "0.14" # for multiplicity dependent efficiency

root4star -b <<EOF
  .L doNbdFitMaker.C
  scan($nevents, "$real", "$mc", $multCut, $nppbin, $nppmin, $nppmax, $kbin, $kmin, $kmax, $xbin, $xmin, $xmax, $efficiency, 1.00, kFALSE);
  .q
EOF

