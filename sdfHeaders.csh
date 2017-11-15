#!/bin/csh
if ( $#argv < 1 ) then
  echo "Converts, fixes headers, and splits out the biases."
  echo "Usage:   prepareIDS.csh [list-of-files]"
  exit
endif
starlink >& /dev/null
convert >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null



foreach file ($argv)
  fitskeys $file:r > zzz_keys.tmp
  set obstype = `grep OBSTYPE zzz_keys.tmp`
  set obstype = $obstype[2]
  echo $obstype  
end
rm zzz_keys.tmp

