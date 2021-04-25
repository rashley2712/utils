#!/bin/csh
if ( $#argv < 1 ) then
  echo "Debiases, and splits the flats and difference science exposures."
  echo "Usage:   prepareFlatsScience.csh  [masterbias]  [list-of-files]"
  exit
endif
starlink >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null

if (! -e bias.reg) then
  cp /storage/astro2/phrnaw/reductions/bias.reg .
  echo "Using the default bias region"
  cat bias.reg
endif 

set bias = $argv[1]
shift

foreach file ($argv)
  $STARLINK_DIR/bin/pamela/debias.pl $bias:r bias.reg $file

  fitskeys $file:r > zzz_keys.tmp
  set obstype = `grep OBSTYPE zzz_keys.tmp`
  set obstype = $obstype[2]
  fitskeys $file:r > zzz_keys.tmp
  set object = `grep "CAT-NAME" zzz_keys.tmp`
  set object = $object[2]

  if ( $obstype == "BIAS" ) then
    if ( ! -e bias ) mkdir bias
    /bin/mv -f $file bias/
    set dir = "bias"
  else if ( $obstype == "FLAT" ) then
    if ( ! -e flat ) mkdir flat
    /bin/mv -f $file flat/
    set dir = "flat"
  else if ( $obstype == "ARC" ) then
    if ( ! -e science ) mkdir science
    /bin/mv -f $file science/
    set dir = "science"
  else
    if ( ! -e science ) mkdir science
    /bin/mv -f $file science
    set dir = "science"
  endif
  echo "File $file has OBSTYPE=$obstype and OBJECT=$object and has been moved to $dir"
end
rm zzz_keys.tmp zzz_debias zzz_debias.log
