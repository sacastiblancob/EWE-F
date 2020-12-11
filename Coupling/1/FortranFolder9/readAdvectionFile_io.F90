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


                 subroutine readAdvectionFile_io &
          (nlat, nlon, Advection_fname)

!  use statevartypesecospace, only: nlat, nlon, grid, Advection_fname, advection
  use statevartypesecospace, only: advection

  implicit none

! INTENT VARIBALES
  integer, INTENT(IN)              :: nlat, nlon
  character(len = 250), INTENT(IN) :: Advection_fname

! IN SUBROUTINE VARIABLES
  integer              :: nrows        ! number of rows
  integer              :: ncols        ! columns of Ecosim input data
  integer              :: iec            ! loop counter
  character ( len = 250 ) :: dummy

  ! open spatial grid data file
  open(400, file = Advection_fname, form = 'formatted')

  ! read number of rows and columns in the file
  read(400, *) dummy, nrows
  read(400, *) dummy, ncols

  ! allocate size of Ecosim data in variable <m_data>
  allocate(advection(nlat, nlon))

  ! read Ecosim input data from file
  do iec = 1, nlat
      read(400, *, end = 100) advection(iec, :)
  end do

100 continue
  close (400)

end subroutine readAdvectionFile_io
