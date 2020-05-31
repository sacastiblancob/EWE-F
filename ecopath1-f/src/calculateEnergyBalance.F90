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

subroutine calculateEnergyBalance ()

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, m_consumed, nvars

implicit none

integer :: var

! calculate production of groups
do var = 1, nvars
    if (ep_data(var)%org_type == 0) then
        ep_data(var)%production = 0
    else
        ep_data(var)%production = ep_data(var)%biomass * ep_data(var)%PoB
    end if
end do

! calculate consumption of consumer groups, else is null
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        ep_data(var)%consumption = ep_data(var)%biomass &
          * ep_data(var)%QoB
    else
        ep_data(var)%consumption = 0
    end if
end do

! predation on groups
do var = 1, nvars
    ep_data(var)%EatenOf = m_consumed(var)
end do

! consumption of groups
do var = 1, nvars
    ep_data(var)%EatenBy = ep_data(var)%consumption
end do

! calculate respiration of consumer groups, else is null
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        ep_data(var)%respiration = ep_data(var)%QoB &
          * ep_data(var)%biomass * (1 - ep_data(var)%unass_Q) &
          - ep_data(var)%PoB * ep_data(var)%biomass
    else
        ep_data(var)%respiration = 0
    end if
end do

! calculate assimilation of consumer groups, else is null
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        ep_data(var)%assimilation = (1 - ep_data(var)%unass_Q) &
          * ep_data(var)%QoB * ep_data(var)%biomass
    else
        ep_data(var)%assimilation = 0
    end if
end do

end subroutine calculateEnergyBalance
