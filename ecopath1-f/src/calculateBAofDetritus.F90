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

subroutine calculateBAofDetritus (ndetritus, detritus_group_no)

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, ep_detfate, m_flows2det, m_consumed

implicit none

integer, intent(in) :: ndetritus, detritus_group_no(ndetritus)
integer             :: i
real(RLEN)          :: Input(ndetritus), surplus


do i = 1, ndetritus
    Input(i) = ep_data(detritus_group_no(i))%detritus_import &
      + sum(m_flows2det(:, i))
end do

do i = 1, ndetritus

    surplus = Input(i) - m_consumed(detritus_group_no(i))

     
    ep_data(detritus_group_no(i))%BA = ep_data(detritus_group_no(i))%BA &
      + surplus * ep_detfate(detritus_group_no(i), i)

end do

end subroutine calculateBAofDetritus
