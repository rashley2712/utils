#!/bin/csh
if ( $#argv < 1 ) then
  echo "Splits BIAS from DATA files. Use after prepareWHT.csh"
  echo "Usage:   splitBiasData.csh [list-of-files]"
  exit
endif
starlink >& /dev/null
convert >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null

if ( ! -e data ) mkdir data
if ( ! -e bias ) mkdir bias

foreach file ($argv)
  fitskeys $file:r > zzz_keys.tmp
  set obstype = `grep OBSTYPE zzz_keys.tmp`
  set obstype = $obstype[2]
  if ( $obstype == "BIAS" ) then
    /bin/mv -f $file bias/
  else
    /bin/mv -f  $file data/
  endif
end

