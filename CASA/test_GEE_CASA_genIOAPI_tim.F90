! these environment variables must be set:
! outname: the I/O API file to write to
! GRIDDESC: description of the grid for I/O API

Program CASA_genIOAPI
  use netcdf
  USE M3UTILIO

  implicit none
  ! include 'netcdf.inc'
  ! include 'PARMS3.EXT'
  ! include 'FDESC3.EXT'
  ! include 'IODECL3.EXT'

  interface
     subroutine get_casa_variable(FILENAME,VARNAME,VAR_DATA,NSTEPS,inlats,inlons)
       character (len=255)::FILENAME! NC PATH/NAME
       character (len=255)::VARNAME! Desired Variable Name
       integer::nsteps
       real,allocatable,intent(inout) :: VAR_DATA(:,:,:)
       real,allocatable,intent(inout) :: inlats(:), inlons(:)
     end subroutine get_casa_variable

     subroutine convert_m2_to_gridcell(in_flux, out_flux, lon, lat)
       REAL, ALLOCATABLE, DIMENSION(:, :, :), INTENT(in) :: in_flux
       REAL, ALLOCATABLE, DIMENSION(:, :, :), INTENT(out) :: out_flux
       REAL, allocatable, INTENT(in) :: lon(:)
       REAL, allocatable, INTENT(in) :: lat(:)
     end subroutine convert_m2_to_gridcell
  end interface

  INTEGER, PARAMETER::&
       NROWS= 181,&
       NCOLS= 288,&
       NSTEPS= 2920,&
       OUT_NROWS = 181,&
       OUT_NCOLS =288,&
       OUT_NSTEPS = 2920,&
       dawn = 13,dusk=3 !UTC
  real, allocatable,dimension(:)::inlats
  real,allocatable,dimension(:)::inlons !one-d arrays of lats and lons
  real, dimension(:, :, :), allocatable::in_co2_flux1, in_co2_flux_per_grid
  real, dimension(out_ncols,out_nrows)::latgrid, longrid !2d arrays of lats and lons
  real, dimension(OUT_NCOLS,OUT_NROWS,nsteps)::CO2_FLUX,casa_temp
  real,dimension(out_ncols,out_nrows)::modis_temp
  real,dimension(3,81600)::modis_ratios
  character(len=255)::FILENAME_modis,varname,filename1
  real::modis_total, d, ratio, finlat, finlon,lat,lon,conversion
  integer::i,j,t,z,i1,casaI,casaJ
  integer::LOGDEV,day,hour,jdate,jtime ,year ,test_loc, ierr ,status

  real::tempd
  real,dimension(out_nsteps)::casasums
  !print*,"making output grids"

  !call make_out_grids(CO2_FLUX,latgrid,longrid,out_nrows,out_ncols,nsteps)
  conversion = -1.0
  print*,"get casa data"
  FILENAME1="/global/project/projectdirs/m2319/Data/CASA/GEE.3hrly.1x1.25.2015.nc"
  VARNAME="GEE"
  ALLOCATE(in_co2_flux1(288, 181, 2920), STAT=ierr)
  ALLOCATE(in_co2_flux_per_grid(288, 181, 2920), STAT=ierr)
  ALLOCATE(inlats(181), STAT=ierr)
  ALLOCATE(inlons(288), STAT=ierr)
  call get_casa_variable(FILENAME1,VARNAME,in_co2_flux1,&
       nsteps,inlats,inlons)

  print*, "casa data loaded!"
  print*, "casa data shape: ", shape(in_co2_flux1)
  print*, "sum CO2(:, :, 1):", sum(in_co2_flux1(:, :, 1))
  print*, "sum CO2(:, :, 2):", sum(in_co2_flux1(:, :, 2))
  print*, "sum CO2(:, :, 2000):", sum(in_co2_flux1(:, :, 2000))
  print*, "write data"
  PRINT*,'STEP4'

  ! Write to ioapi
  nvars3d=1                 ! number of variables
  ftype3d=GRDDED3           ! file is in grided, Global dobson file header
  !gdtyp3d= 2               !Lambert                !
  !gdtyp3d= 1                !lat-lon
  !p_alp3d=33.0              !unused in lat-lon
  !p_bet3d=45.0              !unused in lat-lon
  !p_gam3d=-97.0             !unusued in lat-lon
  !xcent3d=-97.0                !unused in lat-lon
  !ycent3d=39.0                !unused in lat-lon
  !xorig3d=-180
  !yorig3d=-90
  !xcell3d=9000.
  !ycell3d=9000.
  !ncols3d=out_ncols
  !nrows3d=out_nrows
  !nlays3d=1                ! documentation is vague on this maybe vertical levs
  status = DSCGRID( 'CASAGRID', GDNAM3D, GDTYP3D, P_ALP3D, P_BET3D, &
       P_GAM3D, XCENT3D, YCENT3D, XORIG3D, YORIG3D, XCELL3D, YCELL3D, &
       NCOLS3D, NROWS3D, NTHIK3D)
  vgtyp3d=VGSGPN3           !  non-hydrostatic sigma-p vertical coordinate
  vgtop3d=1.                ! domain top in meter
  !nthik3d=1
  NLAYS3D=1
  vglvs3d(1)=1.             ! levels in meter
  units3d(1)='umol/gridcell-1/sec' !the file em_az.dat changes # umol/m2/sec -> molecules/cm2/sec when running STEM
  vname3d(1) = 'CO2_FLUX' !Casa data is in"kgC/m2/s", conversions made on line 110 to convert to umol/m2/sec
  vdesc3d(1) = 'CO2 GEE FLUXES'
  vtype3d(1)=m3real

  sdate3d = 2015001 !also change day below
  stime3d = 0
  tstep3d = 030000
  print*, 'Hello world'
  print*, 'Attempting to open file...'
  if(.not.OPEN3('outname',FSRDWR3,'gen_bdv')) then ! output file does not exit
     print*, 'File does not exist, attempting file creation...'
     if(.not.OPEN3('outname',FSCREA3,'gen_bdv')) then ! FSCREA3 FSUNKN3
        print*, 'Error opening output file'
        stop
     endif
  else
     print*,'Reading from previously created file!'
     if (.not. DESC3('outname') ) then ! if exit, get information
        print*, 'Error getting info from OUTPUT nc'
        stop
     endif
  endif

  ! convert fluxes from per m-2 to per gridcell
  call convert_m2_to_gridcell(in_co2_flux1, in_co2_flux_per_grid, &
       inlons, inlats)
  ! DATA WRITE SECTION
  ! Setup time and iteration information.
  day = 001
  hour = 0
  year = 2015
  jdate =sdate3d
  jtime = stime3d
  ! Attempt an I/O api data write. For whatever reason the convention is to
  !     loop through time and write one 2d grid at a time.
  do t=1,out_nsteps
     ! if(.not.write3('outname',vname3d(1),jdate,jtime,in_co2_flux1(:,:,t))) then
     in_co2_flux_per_grid(:,:,t) = in_co2_flux_per_grid(:,:,t) *83259093
     if(.not.write3('outname',vname3d(1),jdate,jtime,in_co2_flux_per_grid(:,:,t))) then
        print*,'writing error'
        stop
     ENDIF
     print*, 'min in_co2_flux_per_grid this t', minval(in_co2_flux_per_grid(:,:,t))
     call nexttime(jdate,jtime,tstep3d, day, hour,year)
  enddo

  ! do t = 1,out_nsteps
  !    print*, casasums(t), sum(in_co2_flux_per_grid(:,:,t))
  ! enddo
  print*, 'SHUT3()=',SHUT3()
  DEALLOCATE(in_co2_flux1)
  DEALLOCATE(in_co2_flux_per_grid)
  DEALLOCATE(inlats)
  DEALLOCATE(inlons)
end program CASA_genIOAPI

subroutine convert_m2_to_gridcell(in_flux, out_flux, lon, lat)
  USE ioapi_regrid_tools
  implicit NONE

  REAL, allocatable, DIMENSION(:), INTENT(in) :: lon, lat
  REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: area, pct
  REAL, ALLOCATABLE, DIMENSION(:, :, :), INTENT(in) :: in_flux
  REAL, ALLOCATABLE, DIMENSION(:, :, :), INTENT(inout) :: out_flux
  REAL :: cell_EW, cell_NS
  INTEGER status, ierr, dimid, nlon, nlat, ntimes
  INTEGER itime, ilon, ilat, LOGDEV, this_var, this_t, count
  INTEGER jdate, jtime, this_date, this_time

  nlon = SIZE(in_flux, 1)  ! first dimension is longitude
  nlat = SIZE(in_flux, 2)  ! second dimension is latitude
  ntimes = SIZE(in_flux, 3)  ! second dimension is latitude

  !--------------------------------------------------
  ! convert fluxes from umol m-2 s-1 to umol gridcell-1 s-1
  !--------------------------------------------------
  cell_EW = 1.25  ! cell E-W size in degrees longitude

  cell_NS = 1.0   ! cell N-S size in degrees latitude
  ALLOCATE(area(1, 1, nlat, nlon), stat=ierr)
  ALLOCATE(pct(1, 1, nlat, nlon), stat=ierr)

  pct(:, :, :, :) = 1.0  ! no scaling here
  CALL calc_land_area(pct, cell_EW, cell_NS, lon, lat, area)
  ntimes = SIZE(in_flux, 3)  ! third dimension is time
  count = 0
  do itime = 1, ntimes
     do ilon = 1, nlon
        do ilat = 1, nlat
           out_flux(ilon, ilat, itime) = &
                & in_flux(ilon, ilat, itime) * area(1, 1, ilat, ilon)
        ENDDO
     ENDDO
  ENDDO
end subroutine convert_m2_to_gridcell

subroutine nexttime(jdate,jtime,tstep3d,day,hour,year)
  implicit none
  integer::jtime,jdate,tstep3d,day,hour,year
  print*, 'hour 0', hour
  hour = hour+3
  if (hour.eq.24)then
     day = day+1
     hour = 0
  endif
  if (day>365) then
     day = 1
     year = year+1
  endif
  jtime = hour*10000
  jdate = year*1000+day
  return
end subroutine nexttime

subroutine get_casa_variable(FILENAME,VARNAME,VAR_DATA,NSTEPS,inlats,inlons)
  use netcdf
  ! use M3UTILIO
  implicit NONE
  include 'netcdf.inc'
  include 'PARMS3.EXT'
  include 'FDESC3.EXT'
  include 'IODECL3.EXT'

  REAL, PARAMETER::&
       XCELL = 1.25,&
       YCELL = 1.0
  integer::nsteps, this_t, ierr
  character (len=255)::FILENAME! NC PATH/NAME
  character (len=255)::VARNAME! Desired Variable Name
 ! character (len=*),parameter::&
 !      lat_name="latitude",&
 !      lon_name="longitude"
  integer::i,j,t! indexes for iterating
  integer::ncid,varid!id for nc reading
  integer::numDims,numAtts,RecordDimLength
  character(len = nf90_max_name) :: RecordDimName
  integer, dimension(nf90_max_var_dims)::DimIds
  REAL, allocatable, DIMENSION(:, :, :), intent(inout) :: var_data
  ! real, allocatable, dimension(288,181,nsteps)::var_data !data array
  real, dimension(288, 181, 1)::var_data_onestep !data array

  real,allocatable, dimension(:),intent(inout)::inlats
  real,allocatable,dimension(:),intent(inout)::inlons
  !real::templat,templon

  !       print*,"enter get_casa_variables()"
  print*,"Open: ",FILENAME
  call check(nf90_open(FILENAME,NF90_NOWRITE,ncid))

  ! Read actual data
  print*,"Find variable: ",VARNAME
  call check(nf90_inq_varid(ncid,VARNAME,varid))
  print*,"Load variable CASA: ",VARNAME
  print*,'var_data shape: ', shape(var_data), 'nsteps: ', nsteps
  call check(nf90_inquire_variable(ncid, varid, ndims = numDims, &
       natts = numAtts, dimids = DimIds))
  print*,'nf90_inquire_variable: ndims:', numDims, 'natts', numAtts
  call check(nf90_inquire_dimension(ncid, 1, RecordDimName, RecordDimLength))
  print*,'nf90_inquire_dimension 1:', RecordDimName, RecordDimLength
  call check(nf90_inquire_dimension(ncid, 2, RecordDimName, RecordDimLength))
  print*,'nf90_inquire_dimension 2:',  RecordDimName, RecordDimLength
  call check(nf90_inquire_dimension(ncid, 3, RecordDimName, RecordDimLength))
  print*,'nf90_inquire_dimension 3:', RecordDimName, RecordDimLength
  call check(nf90_get_var(ncid,varid,var_data))
  print*, "sum GPP CO2(:, :, 1) (in subroutine):", sum(var_data(:, :, 1))
  print*, "sum GPP CO2(:, :, 2) (in subroutine):", sum(var_data(:, :, 2))
  print*, "sum GPP CO2(:, :, 2000) (in subroutine):", sum(var_data(:, :, 2000))
  print*, "load inlats and inlons"
  call check(nf90_inq_varid(ncid,"lon",varid))
  call check(nf90_get_var(ncid,varid,inlons))
  print*, "Longitude loaded!"
  call check(nf90_inq_varid(ncid,"lat",varid))
  call check(nf90_get_var(ncid,varid,inlats))
  print*, "Latitude Loated!"

  print*,"range of CASA GPP co2",maxval(var_data),minval(var_data)
  call check(nf90_close(ncid))
  print*,"range of CASA GPP co2",maxval(var_data),minval(var_data)
  print*,"File closed."
  ! load longitude and latitude arrays

  !reshape the data from the buffer into VAR_DATA
  return
end subroutine get_casa_variable

subroutine check(status)
  use netcdf
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
end subroutine check
