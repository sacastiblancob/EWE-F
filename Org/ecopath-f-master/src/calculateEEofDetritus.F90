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

subroutine calculateEEofDetritus (ndetritus, detritus_group_no)

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, m_flows2det, m_consumed

implicit none

integer, intent(in) :: ndetritus, detritus_group_no(ndetritus)
integer             :: i


! calculate EE of detritus groups
do i = 1, ndetritus
    ep_data(detritus_group_no(i))%EE = m_consumed(detritus_group_no(i)) &
      / (sum(m_flows2det(:, i)) &
      + ep_data(detritus_group_no(i))%detritus_import)
end do

end subroutine calculateEEofDetritus
