#!/bin/csh
if ( $#argv < 1 ) then
  echo "Usage: makebias [list-of-files]"
  exit
endif
starlink >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null

if ( -e bias.lis ) rm bias.lis
foreach file ($argv)
  echo $file:r >> bias.lis
end
medsky bias.lis masterbias scaled=FALSE
gaia masterbias

