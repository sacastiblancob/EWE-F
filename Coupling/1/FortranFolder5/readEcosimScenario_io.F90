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

            subroutine readEcosimScenario_io &
      (GroupInfo_fname, nvars, nstanzas)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

! This subroutine reads Group Info paramters of Ecosim
! from tab-delimited ascii file

use statevartypesecosim, only: groupnames, es_data
!use statevartypesecosim, only: GroupInfo_fname, nvars, nstanzas, groupnames

implicit none

!INTENT VARIABLES
character(len = 250), intent(in)                 :: GroupInfo_fname
integer, intent(inout)                           :: nvars, nstanzas
!character(len = 250), allocatable, intent(inout) :: groupnames(:)
!intent(inout)                                    :: es_data

!IN SUBROUTINE VARIABLES
!integer                :: nrows        ! number of rows
integer                 :: ncols        ! columns of Ecosim input data
integer                 :: iec            ! loop counter
real(RLEN), allocatable :: m_data(:, :) ! matrix to write into Ecosim input
character ( len = 250 ) :: dummy


! open Ecosim data file
open(111, file = GroupInfo_fname, form = 'formatted')

! read number of rows and columns in the file
read(111, *) dummy, nvars
read(111, *) dummy, ncols
read(111, *) dummy, nstanzas
read(111, *)

! number of state variables <nvars> is equal to the number of rows
! <nrows> in the scenario parameters matrix
!nvars = nrows

! allocate size of Ecosim data in variable <m_data>
allocate(m_data(nvars, ncols))
allocate(groupnames(nvars))

! read Ecosim input data from file
do iec = 1, nvars
    read(111, *, end = 100) groupnames(iec), m_data(iec, :)
end do

100 continue
 close (111)

! allocate and populate fields of <es_data>
! call setEcosimScenarioParameters(nrows, ncols, m_data, es_data, nvars)

 allocate(es_data(nvars))

do iec = 1, nvars
    ! populate fields of es_data(iec) from scenario parameters
    es_data(iec)%rel_PoB_max  = m_data(iec, 1)
    es_data(iec)%Ftime        = m_data(iec, 2)
    es_data(iec)%Ftime_max    = m_data(iec, 3)
    es_data(iec)%Ftime_adjust = m_data(iec, 4)
    es_data(iec)%M0_pred      = m_data(iec, 5)
    es_data(iec)%risk_time    = m_data(iec, 6)
    es_data(iec)%Q_maxoQ_0    = m_data(iec, 7)
    es_data(iec)%QB_maxoQB_0  = m_data(iec, 8)
    es_data(iec)%switch_power = m_data(iec, 9)
    es_data(iec)%isAdvected   = m_data(iec, 10)
end do

 deallocate(m_data)

end subroutine readEcosimScenario_io
