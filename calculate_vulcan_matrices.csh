#!/bin/csh

####SBATCH -p debug
####SBATCH -N 1
####SBATCH -t 00:01:00
####SBATCH -L project     #note: specify license need for the file systems your job needs, such as SCRATCH,project
####SBATCH -o VULCAN.%A.out

setenv SCRATCHDIR /global/cscratch1/sd/twhilton/Vulcan
setenv COL_REFINEMENT  3
setenv ROW_REFINEMENT  3
setenv SCALEFAC        1.0
setenv GRIDDESC        $HOME/Code/Regrid/VulcanRegrid/GRIDDESC_VULCAN
setenv MATRIX          $SCRATCHDIR/VULCANto9kmSTEMmatrix
setenv MATTXT          $SCRATCHDIR/VULCANto9kmSTEMmattxt

rm -fv $MATRIX
rm -fv $MATTXT

setenv PROMPTFLAG TRUE
echo "starting mtxcalc: "
date
mtxcalc << DONE
Y
VULCANGRID
STEM_9KM_GRD
MATTXT
DONE
echo "finished mtxcalc: "
date

setenv IN_DATA $SCRATCHDIR/reversed_vulcan_ioapi_1jan.nc
setenv OUT_DATA $SCRATCHDIR/vulcan_fossilCO2_stem9km_ioapi.nc

rm -fv $OUT_DATA

echo "starting mtxcple: "
date
mtxcple << DONE
Y
NONE
MATRIX
IN_DATA




1
OUT_DATA
DONE
echo "finished mtxcalc: "
date
