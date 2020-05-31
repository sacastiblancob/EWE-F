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

subroutine setEcosimScenarioParameters (nrows, ncols, m_data)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

! This subroutine populates the fields of input data <es_data>
! after scenario parameters and initial conditions are
! read from the file into temporary
! matrix named <m_data>. The type of the columns of the m_data
! matrix are defined in statevartypesecosim.f90 and is fixed in order.
! However, they can be extended by appending new fields to <es_data>.

  use statevartypesecosim, only: es_data, nvars

implicit none

! variables inherited from ecopath model
integer, intent(in) :: nrows, ncols
real(RLEN), intent(in) :: m_data(nrows, ncols)

! in-subroutine variables
integer             :: i      ! loop variable

allocate(es_data(nvars))

do i = 1, nrows
    ! populate fields of es_data(i) from scenario parameters
    es_data(i)%rel_PoB_max  = m_data(i, 1)
    es_data(i)%Ftime        = m_data(i, 2)
    es_data(i)%Ftime_max    = m_data(i, 3)
    es_data(i)%Ftime_adjust = m_data(i, 4)
    es_data(i)%M0_pred      = m_data(i, 5)
    es_data(i)%risk_time    = m_data(i, 6)
    es_data(i)%Q_maxoQ_0    = m_data(i, 7)
    es_data(i)%QB_maxoQB_0  = m_data(i, 8)
    es_data(i)%switch_power = m_data(i, 9)
    es_data(i)%isAdvected   = m_data(i, 10)
end do

end subroutine setEcosimScenarioParameters
