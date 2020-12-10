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

subroutine readEcopathScenario_io ()

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: InputData_fname, isASCIIinputFile, &
                          nvars, nstanzas, groupnames

implicit none

integer              :: nrows        ! number of rows
integer              :: ncols        ! columns of Ecopath input data
integer              :: i            ! loop counter
real(RLEN), allocatable :: m_data(:, :) ! matrix to write into Ecopath input
character ( len = 250 ) :: dummy

! open Ecopath data file
open(333, file = InputData_fname, form = 'formatted' )

! read number of rows and columns in the file
read(333, *) dummy, nrows
read(333, *) dummy, ncols
read(333, *) dummy, nstanzas
read(333, *)

! number of state variables <nvars> is equal to the number of rows 
! <nrows> in the scenario parameters matrix
nvars = nrows

! allocate size of Ecopath data in variable <m_data>
allocate(m_data(nrows, ncols))
allocate(groupnames(nrows))

! read Ecopath input data from file
if (isASCIIinputFile) then
    do i = 1, nrows
        read(333, *, end = 300) groupnames(i), m_data(i, :)
    end do
end if

300 continue
 close(333)

! allocate and populate fields of <ep_data>
if (isASCIIinputFile) then
    call setEcopathScenarioParameters (nrows, ncols, m_data)
end if

deallocate(m_data)

end subroutine readEcopathScenario_io
