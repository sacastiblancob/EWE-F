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

subroutine calculateDetritusFate (ndetritus, detritus_group_no)

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, ep_detfate, m_flows2det, &
                          m_consumed

implicit none

integer, intent(in) :: ndetritus, detritus_group_no(ndetritus)
integer             :: i, j
real(RLEN)             :: Input(ndetritus), DetPassedOn(ndetritus), surplus

ep_data(:)%DetPassedProp = 0

do i = 1, ndetritus
    Input(i) = ep_data(detritus_group_no(i))%detritus_import &
      + sum(m_flows2det(:, i))
end do

do i = 1, ndetritus

    DetPassedOn(i) = 0
    surplus        = Input(i) - m_consumed(detritus_group_no(i))

    if (surplus > 0) then
        do j = 1, ndetritus
            if (i /= j) then
                Input(j) = Input(j) + surplus &
                  * ep_detfate(detritus_group_no(i), j)

                DetPassedOn(i) = DetPassedOn(i) + surplus &
                  * ep_detfate(detritus_group_no(i), j)

                m_flows2det(detritus_group_no(i), j) = surplus &
                  * ep_detfate(detritus_group_no(i), j)
            end if
        end do
    end if

    ep_data(detritus_group_no(i))%DetPassedProp = real(DetPassedOn(i) &
      / (ep_data(detritus_group_no(i))%biomass + 1.0e-20), 4)

end do

end subroutine calculateDetritusFate
