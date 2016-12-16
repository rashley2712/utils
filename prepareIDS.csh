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


if ( ! -e fitsfiles ) mkdir fitsfiles
if ( ! -e data ) mkdir data
if ( ! -e bias ) mkdir bias

foreach file ($argv)
  chmod +w $file
  set sdffile = $file:r".sdf"
  if ( ! -e $sdffile ) then
    fits2ndf $file $file:r
    echo "$file  --> $sdffile"
  else
    echo "$sdffile already exists"
  endif
  mv $file fitsfiles
  ultradas $sdffile
  fixhead 1 7 $sdffile
  fitskeys $sdffile:r > zzz_keys.tmp
  set obstype = `grep OBSTYPE zzz_keys.tmp`
  set obstype = $obstype[2]
  if ( $obstype == "BIAS" ) then
    /bin/mv -f $sdffile bias/
  else
    /bin/mv -f  $sdffile data/
  endif
end
rm zzz_keys.tmp zzz_fixhead zzz_fixhead.log

