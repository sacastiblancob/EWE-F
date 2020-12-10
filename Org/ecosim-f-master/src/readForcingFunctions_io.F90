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

subroutine readForcingFunctions_io ()

! This subroutine reads forcing time series data
! from tab-delimited ascii file.

use statevartypesecosim, only: Forcing_fname, force, nvars

implicit none

integer              :: i             ! loop counter
integer              :: frows, fcols  ! number of rows and columns 
integer, allocatable :: groupno(:)    ! number of groups in forcing data
character ( len = 250 ) :: groupname, dummy

! open time-series data file
open(555, file = Forcing_fname, form = 'formatted')

allocate(groupno(nvars))
allocate(force%forcetype(nvars))

! read number of group, time-series type, 
! number of rows and columns in file
read(555, *) dummy, frows
read(555, *) dummy, fcols
read(555, *) dummy, groupno(:)
read(555, *) dummy, force%forcetype(:)
read(555, *) groupname

! allocate the time-series to their corresponding groups
allocate(force%fishforce(frows, fcols))

! read time-series data from file
do i = 1, frows
    read(555, *, end = 500) dummy, force%fishforce(i, :)
end do

500 continue
 close(555)

end subroutine readForcingFunctions_io
