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

subroutine calculateDetritalFlows (detritus_no)

! This subroutine calculate flows to detritus compartments
! only from living groups. Flows between multiple detritus groups
! are handled elsewhere.

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, ep_detfate, m_flows2det, nvars

implicit none

integer             :: var           ! loop counters
integer, intent(in) :: detritus_no   ! detritus number

! calculate the flows from living groups to each detritus compartment
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        ! group is consumer, account for unassimilated consumption
        m_flows2det(var, detritus_no) = (ep_data(var)%biomass &
          * ep_data(var)%PoB * (1 - ep_data(var)%EE) &
          + (ep_data(var)%biomass * ep_data(var)%QoB &
          * ep_data(var)%unass_Q)) * ep_detfate(var, detritus_no) &
          + ep_data(var)%discards

    else if (ep_data(var)%org_type == 1) then
        ! group is producer than there is only natural mortality
        m_flows2det(var, detritus_no) = (ep_data(var)%biomass &
          * ep_data(var)%PoB * (1 - ep_data(var)%EE)) &
          * ep_detfate(var, detritus_no) &
          + ep_data(var)%discards
     else
        ! this group is detritus
        continue          
    end if
end do

end subroutine calculateDetritalFlows
