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

              subroutine readSpatialDistribution_io &
       (nlat, nlon, nvars, SpatialDistribution_dirname)

!  use statevartypesecosim, only: GroupInfo_fname, nvars
  use statevartypesecospace, only: spatialhafs

  implicit none

! INTENT VARIABLES
  INTEGER, INTENT(IN)               :: nlat, nlon, nvars
  character(len = 250), INTENT(IN)  :: SpatialDistribution_dirname

! IN SUBROUTINE VARIABLES
  integer :: iec, j
!  integer              :: nrows        ! number of rows
  character(len = 250) :: fileinptr
!  character ( len = 250 ) :: dummy

! added by Sergio to solve hardcoded 22 for nvars
! open Ecosim data file
!  open(6666, file = GroupInfo_fname, form = 'formatted')

! read number of rows
!  read(6666, *) dummy, nrows
!  close(6666)

! number of state variables <nvars> is equal to the number of rows
! <nrows> in the scenario parameters matrix
!  nvars = nrows

! nvars = 22;	!original

  allocate(spatialhafs(nlat, nlon, nvars))

  do iec = 1, nvars
      if (iec < 10) then
          write(fileinptr, '(A, I1, ".csv")') TRIM(SpatialDistribution_dirname), iec
          print *,fileinptr
      else
          write(fileinptr, '(A, I2, ".csv")') TRIM(SpatialDistribution_dirname), iec
          print *,fileinptr
      end if

      print *, fileinptr
      open(4440444, file = fileinptr, form = 'formatted')
      read(4440444, *)
      read(4440444, *)

      do j = 1, nlat
          read(4440444, *, end = 1000) spatialhafs(j, :, iec)
      end do

1000  continue
  close(4440444)

  end do

end subroutine readSpatialDistribution_io
