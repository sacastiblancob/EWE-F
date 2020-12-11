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

subroutine calculateMaximumPoBRelatedValues(nvars, ep_data, es_data)

use statevartypesecopath, only: ecopath_data
use statevartypesecosim, only: ecosim_data

implicit none

! INTENT VARIABLES
INTEGER, INTENT(IN)              :: nvars
TYPE(ecopath_data), INTENT(IN)   :: ep_data(nvars)
TYPE(ecosim_data), INTENT(INOUT) :: es_data(nvars)

! IN SUBROUTINE VARIABLES
integer :: var


do var = 1, nvars
    if (ep_data(var)%org_type == 0 .or. ep_data(var)%org_type == 2) then
        es_data(var)%abs_PoB_max = 0
    else
        es_data(var)%abs_PoB_max = es_data(var)%rel_PoB_max &
          * ep_data(var)%PoB
!        write(*,*) es_data(var)%abs_PoB_max
    end if
end do

do var = 1, nvars
    if (ep_data(var)%org_type == 1) then
        es_data(var)%PoB_biomass = (es_data(var)%abs_PoB_max &
          / ep_data(var)%PoB - 1) / ep_data(var)%biomass
!        write(*,*) es_data(var)%PoB_biomass
    else
        es_data(var)%PoB_biomass = 0
    end if
!    write(*,*) ep_data(var)%biomass
end do

end subroutine calculateMaximumPoBRelatedValues
