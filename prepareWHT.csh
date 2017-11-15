#!/bin/csh
if ( $#argv < 1 ) then
  echo "Fixes headers and sorts out the blue-arm and red-arm data"
  echo "Usage: prepareWHT [list-of-files]"
  exit
endif
starlink >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null
convert >& /dev/null

if ( ! -e red  ) mkdir red
if ( ! -e blue ) mkdir blue
if ( ! -e aux  ) mkdir aux
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
  fixhead 1 3 $sdffile
  fitskeys $sdffile:r > zzz_keys.tmp
  set detector = `grep DETECTOR zzz_keys.tmp`
  set detector = $detector[2]
  if ( $detector == "EEV12" ) mv $sdffile:r".sdf" blue/
  if ( $detector == "REDPLUS" ) mv $sdffile:r".sdf" red/
  if ( $detector == "AUXCAM" ) mv $sdffile:r".sdf" aux/
end
rm zzz_fixhead* zzz_keys.tmp

