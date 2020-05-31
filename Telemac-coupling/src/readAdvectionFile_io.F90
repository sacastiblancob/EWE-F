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

 
subroutine readAdvectionFile_io ()

  use statevartypesecospace, only: nlat, nlon, grid, Advection_fname, advection
  
  implicit none

  integer              :: nrows        ! number of rows
  integer              :: ncols        ! columns of Ecosim input data
  integer              :: i            ! loop counter
  character ( len = 250 ) :: dummy

  ! open spatial grid data file
  open(400, file = Advection_fname, form = 'formatted')

  ! read number of rows and columns in the file
  read(400, *) dummy, nrows
  read(400, *) dummy, ncols

  nlat = nrows
  nlon = ncols

  ! allocate size of Ecosim data in variable <m_data>
  allocate(advection(nlat, nlon))

  ! read Ecosim input data from file
  do i = 1, nlat
      read(400, *, end = 100) advection(i, :)
  end do

100 continue
  close (400)

end subroutine readAdvectionFile_io
