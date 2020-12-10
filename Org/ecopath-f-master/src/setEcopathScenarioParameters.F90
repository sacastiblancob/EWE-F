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

subroutine setEcopathScenarioParameters (nrows, ncols, m_data)

! This subroutine populates the fields of input data <ep_data>
! after input parameters are read from the file into a temporary
! matrix named <m_data>. The type of the columns of the m_data
! matrix are defined in statevartypesecopath.f90.
! However, they can be extended by appending new fields to <ep_data>.

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, nvars

implicit none

! variables inherited from Ecopath-FORTRAN model
integer, intent(in) :: nrows, ncols
real(RLEN), intent(in) :: m_data(nrows, ncols)
integer             :: i

! allocate and populate fields of ep_data
allocate(ep_data(nvars))

do i = 1, nrows
    ep_data(i)%biomass = m_data(i, 1)
    ep_data(i)%PoB = m_data(i, 2)
    ep_data(i)%QoB = m_data(i, 3)
    ep_data(i)%EE = m_data(i, 4)
    ep_data(i)%PoQ = m_data(i, 5)
    ep_data(i)%unass_Q = m_data(i, 6)
    ep_data(i)%detritus_import = m_data(i, 7)
    ep_data(i)%landings = m_data(i, 8)
    ep_data(i)%discards = m_data(i, 9)
    ep_data(i)%BA = m_data(i, 10)
    ep_data(i)%org_type = int(m_data(i, 11))
    ep_data(i)%isstanza = int(m_data(i, 12))
    ep_data(i)%stanza_no = int(m_data(i, 13))
    ep_data(i)%age_start = int(m_data(i, 14))
    ep_data(i)%isleading = int(m_data(i, 15))
end do

end subroutine setEcopathScenarioParameters
