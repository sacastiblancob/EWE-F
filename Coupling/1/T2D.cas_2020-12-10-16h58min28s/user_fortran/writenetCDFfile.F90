!========================================================================
! This file is part of EwE-F
! Copyright (C) 2011-2019 Middle East Technical University
! Institute of Marine Sciences (IMS-METU), Erdemli/Turkey and
! Istituto Nazionale di Oceanografia e di Geofisica Sperimentale (OGS),
! Trieste/Italy.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 as
! published by the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with  this program; if not,
! see <http://www.gnu.org/licenses/gpl-2.0.html>.
!========================================================================

subroutine writenetCDFfile (noftsteps, mat_out, ncdfout_fname)
  use netcdf
  use statevartypesecosim, only: nvars, RLEN
  use statevartypesecospace, only: nlat, nlon
  implicit none

  ! This is the name of the data file we will create.
  character (len = *), intent(in)          :: ncdfout_fname
  integer                        :: ncid
  integer, intent(in)            :: noftsteps

  ! We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
  ! timesteps of data.
  integer, parameter             :: NDIMS = 4
  character (len = *), parameter :: LVL_NAME = "functionalGroup"
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: REC_NAME = "time"
  integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid
  integer :: NRECS, NLVLS, NLATS, NLONS


  ! These program variables hold the latitudes and longitudes.
  integer :: mat_varid

  ! We will create two netCDF variables, one each for temperature and
  ! pressure fields.
  integer :: dimids(NDIMS)

  ! Program variables to hold the data we will write out. We will only
  ! need enough space to hold one timestep of data; one record.
  real(RLEN), intent(in) :: mat_out(nlat, nlon, noftsteps + 1, nvars)

  NRECS = noftsteps + 1
  NLVLS = nvars
  NLATS = nlat
  NLONS = nlon

  WRITE(*,*) ncdfout_fname

  ! Create the file.
  call check( nf90_create(ncdfout_fname, nf90_clobber, ncid) )

  WRITE(*,*) ncdfout_fname

  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
  call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
  call check( nf90_def_dim(ncid, REC_NAME, NRECS, rec_dimid) )
  call check( nf90_def_dim(ncid, LVL_NAME, NLVLS, lvl_dimid) )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids = (/ 1,2,3,4 /)

  ! Define the variable
  call check(nf90_def_var(ncid, "Biomass", NF90_FLOAT, dimids, mat_varid))
  call check(nf90_put_att(ncid, mat_varid, "_FillValue", real(-999)))

  ! End define mode.
  call check( nf90_enddef(ncid) )

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays only hold one timestep worth
  ! of data. We will just rewrite the same data for each timestep. In
  ! a real :: application, the data would change between timesteps.


  call check( nf90_put_var(ncid, mat_varid, mat_out) )

  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are really written to disk.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS writing example file ", ncdfout_fname, "!"

contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end subroutine writenetCDFfile
