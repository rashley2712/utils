#!/bin/csh
if ( $#argv < 1 ) then
  echo "Show the headers for a bunch of SDF files"
  echo "Usage:   showheader.csh [header] [list-of-files]"
  exit
endif
starlink >& /dev/null
convert >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null

set header = $argv[1]
shift

foreach file ($argv)
  fitskeys $file:r > zzz_keys.tmp
  set value = `grep $header zzz_keys.tmp`
  set value = $value[2]
  echo $file $header $value
  # cat zzz_keys.tmp | grep $header
end

