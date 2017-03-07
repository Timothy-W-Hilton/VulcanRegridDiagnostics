program CASA_gee_per_m2
  ! create an I/O API file containing CASA GEE in umol C m-2 s-1 from
  ! an I/O API file containing CASA GEE in umol C cell-1 s-1.
  !
  ! to compile: use makefile casa_2_m2.mk:
  !     make -f casa_2_m2.mk clobber all
  !
  ! Timothy W. Hilton, thilton@ucmerced.edu, 7 Mar 2017
  USE ioapi_regrid_tools

  IMPLICIT NONE

  character*80 vname_arg
  character*80 units_arg
  character*80 vdesc_arg
  real m2_per_cell

  vname_arg = 'CO2_FLUX'
  units_arg = 'umol C m-2 s-1'
  vdesc_arg = 'CASA GEE'

  m2_per_cell = 9000 * 9000 ! gridcell size in m2: 9 km2 per gridcell

  ! 1. need to set INPUT and OUTPUT logical names for this call
  ! 2. per_cell_2_per_m2 is defined in ioapi_regrid_tools.F90, which
  ! is compiled to the library
  ! /project/projectdirs/m2319/local/lib/libioapi_regrid_tools.a
  CALL per_cell_2_per_m2(vname_arg, units_arg, vdesc_arg, m2_per_cell)

END program CASA_gee_per_m2
