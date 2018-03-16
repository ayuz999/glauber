#!/bin/csh -f

#set energy = "11"
set energy = "22"

#set list = ( "default" "small" "large" "smallXsec" "largeXsec" "gray" "gauss" )
set list = ( "default" "small" "large" "smallXsec" "largeXsec" "gauss" "smallNpp" "largeNpp" )
#set list =("default")


foreach type ($list)
  set output = "LIST/tree.$type.list"
  echo "Make $output ..."
#  ls -1 outputs/icmaker_Au_39_${type}_spherical*root > $output
   ls -1 output/fastglaubermc_CuCu_${energy}*GeV_${type}_spherical_*.root > $output
end

# Make link "smallNpp" "largeNpp" "smallTotal" "largeTotal" to "default"
cd LIST
ln -vs tree.default.list tree.smallNpp.list
ln -vs tree.default.list tree.largeNpp.list
ln -vs tree.default.list tree.smallTotal.list
ln -vs tree.default.list tree.largeTotal.list
ln -vs tree.default.list tree.lowrw.list
ln -vs tree.default.list tree.highrw.list
ls -lrta
cd ../

