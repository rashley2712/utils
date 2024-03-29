#!/bin/csh
if ( $#argv < 2 ) then
  echo "Usage: makeflat.csh [polynomial-order] [mask? y/n] [list-of-files]"
  exit
endif
starlink >& /dev/null
kappa >& /dev/null
figaro >& /dev/null
pamela >& /dev/null
if ( -e flat.lis ) rm flat.lis
set npoly = $argv[1]
shift
set mask = $argv[1]
shift

foreach file (`cat $argv`)
  echo $file:r >> flat.lis
end

medsky flat.lis mflat scaled=TRUE

ystract mflat min max mflaty
log10 mflaty mflatylog

if ( $mask == "y" ) then
  echo "Mask mask.ard will be applied"
  ardmask mflatylog mask.ard mflatylogmask
  mv -f mflatylogmask.sdf mflatylog.sdf
endif

polfit mflatylog mflatylogfit $npoly -3 3 50
exp10 mflatylogfit mflatyfit
isydiv mflat mflatyfit mflatfit

istat mflatfit 1600 2000 70 250 > zzz_istat.tmp
set mean = `grep Mean zzz_istat.tmp`
set mean = $mean[3]

icdiv mflatfit $mean mflatnorm
isub mflatnorm mflatnorm zero
icadd zero 1 unit
idiv unit mflatnorm balance

rm zzz_istat.tmp
gaia balance.sdf

