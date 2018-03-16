#!/bin/csh
set nppbin  = "1"
set nppmin  = "2.20"
set nppmax  = "2.20"

set xbin    = "1"
set xmin    = "0.13"
set xmax    = "0.13"

if ( $#argv != 1 ) then
  echo ""
  echo " Usage : $0 [Input k value]"
  echo ""
  echo "    Scan data for npp = $nppmin - $nppmax with input x = $xmin - $xmax and k value"
  echo ""
  exit
endif

starver SL16a

set kbin    = "1"
set kmin    = "$1"
set kmax    = "$1"

set nevents  = "1000000" # 1M
set real     = "refmult3.root"
#set real     = "npart.root"
set mc       = "ncoll_npart.root"
set multCut  = "100" #for low multiplicity cut
set efficiency = "0.14" # for multiplicity dependent efficiency

root4star -b <<EOF
  .L doNbdFitMaker.C
  scan($nevents, "$real", "$mc", $multCut, $nppbin, $nppmin, $nppmax, $kbin, $kmin, $kmax, $xbin, $xmin, $xmax, $efficiency, 1.00, kFALSE);
  .q
EOF

