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


               subroutine readGridFile_io &
          (nlat, nlon, SpatialGrid_fname)

  use statevartypesecospace, only: grid
!  use statevartypesecospace, only: nlat, nlon, grid, SpatialGrid_fname

  implicit none

! INTENT VARIABLES
  INTEGER, INTENT(INOUT)                 :: nlat, nlon
  character(len = 250), INTENT(IN)       :: SpatialGrid_fname
!  integer, allocatable, INTENT(OUT)    :: grid(:, :)

! IN SUBROUTINE VARIABLES
!  integer              :: nrows        ! number of rows
!  integer              :: ncols        ! columns of Ecosim input data
  integer              :: i            ! loop counter
  character ( len = 250 ) :: dummy

  ! open spatial grid data file
  open(400, file = SpatialGrid_fname, form = 'formatted')

  ! read number of rows and columns in the file
  read(400, *) dummy, nlat
  read(400, *) dummy, nlon

!  nlat = nrows
!  nlon = ncols

  ! allocate size of Ecosim data in variable <m_data>
  allocate(grid(nlat, nlon))

  ! read Ecosim input data from file
  do i = 1, nlat
      read(400, *, end = 100) grid(i, :)
  end do

100 continue
  close (400)

end subroutine readGridFile_io
