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

subroutine readEcosimScenario_io ()

#ifdef isWithBFM
  use global_mem
#else
  use statevartypesecopath, only: RLEN
#endif

! This subroutine reads Group Info paramters of Ecosim
! from tab-delimited ascii file

use statevartypesecosim, only: GroupInfo_fname, nvars, nstanzas, groupnames

implicit none

integer              :: nrows        ! number of rows
integer              :: ncols        ! columns of Ecosim input data
integer              :: i            ! loop counter
real(RLEN), allocatable :: m_data(:, :) ! matrix to write into Ecosim input
character ( len = 250 ) :: dummy

! open Ecosim data file
open(111, file = GroupInfo_fname, form = 'formatted')

! read number of rows and columns in the file
read(111, *) dummy, nrows
read(111, *) dummy, ncols
read(111, *) dummy, nstanzas
read(111, *)

! number of state variables <nvars> is equal to the number of rows 
! <nrows> in the scenario parameters matrix
nvars = nrows

! allocate size of Ecosim data in variable <m_data>
allocate(m_data(nrows, ncols))
allocate(groupnames(nrows))

! read Ecosim input data from file
do i = 1, nrows
    read(111, *, end = 100) groupnames(i), m_data(i, :)
end do

100 continue
 close (111)

! allocate and populate fields of <es_data>
 call setEcosimScenarioParameters (nrows, ncols, m_data)

 deallocate(m_data)

end subroutine readEcosimScenario_io
