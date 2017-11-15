#!/bin/csh
if ( $#argv < 1 ) then
  echo "Removes acquisition images from the folder."
  echo "Usage:   delacquisition.csh [list-of-files]"
  exit
endif
starlink >& /dev/null
convert >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null

if ( ! -e acquisition ) mkdir acquisition

foreach file ($argv)
  fitskeys $file:r > zzz_keys.tmp
  set obstype = `grep OBSTYPE zzz_keys.tmp`
  echo OBSTYPE: $obstype
  if ("$obstype" == "") then
    echo "OBSTYPE is blank"
    /bin/mv -f $file acquisition/

  endif
end

