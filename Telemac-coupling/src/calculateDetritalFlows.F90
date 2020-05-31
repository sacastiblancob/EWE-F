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

subroutine calculateDetritalFlows (biomass)

!#ifdef isWithBFM
!  use global_mem
!  use mem, ONLY: iiHigherTrophicLevels
!#else
  use statevartypesecopath, only: RLEN
!#endif

use statevartypesecopath, only: ep_data, ep_detfate, flow2detritus, &
                                det_export_rate
use statevartypesecosim, only: es_data, nvars, FirstTime, detritus_no, &
                               ndetritus

implicit none

! variables inherited from Ecosim-FORTRAN model
real(RLEN), intent(in) :: biomass(nvars)
!#ifdef isWithBFM
!real(RLEN)             :: reR6c(1)
!#endif

! in-subroutine variables
integer :: i, j, var

!!!!! calculate detritus export rate

do j = 1, ndetritus

    flow2detritus(j) = 0

!#ifdef isWithBFM
!    do var = 1, iiHigherTrophicLevels
!        reR6c = 0.
!#endif
!#ifndef isWithBFM
    do var = 1, nvars
!#endif
        if (ep_data(var)%org_type /= 0) then
            ! calculate non-predation natural mortality
            ! for groups other than detritus
!#ifdef isWithBFM
!            reR6c(1) = ((1 - ep_data(var)%EE) &
!              * ep_data(var)%PoB) * biomass(var) * ep_detfate(var, j)
!#endif
            flow2detritus(j) = flow2detritus(j) + ((1 - ep_data(var)%EE) &
              * ep_data(var)%PoB) * biomass(var) * ep_detfate(var, j)
        end if

        if (ep_data(var)%org_type == 2) then
            ! if group is a consumer, consider calculating excretion
!#ifdef isWithBFM
!            reR6c(1) = reR6c(1) + ep_data(var)%unass_Q &
!              * es_data(var)%EatenBy * ep_detfate(var, j)
!#endif
            flow2detritus(j) = flow2detritus(j) + ep_data(var)%unass_Q &
              * es_data(var)%EatenBy * ep_detfate(var, j)
        end if

        if (ep_data(var)%discards /= 0) then
            if (FirstTime .eqv. .true.) then
                flow2detritus(j)  = flow2detritus(j) + (ep_data(var)%landings &
               + ep_data(var)%discards) * (ep_data(var)%discards &
               / (ep_data(var)%landings + ep_data(var)%discards))
            else
                ! then calculate flows from discarded fish
                flow2detritus(j)  = flow2detritus(j) + biomass(var) &
                  * es_data(var)%fishmort * (ep_data(var)%discards &
                  / (ep_data(var)%landings + ep_data(var)%discards))
            end if
        end if

!#ifdef isWithBFM
!    call calculateMultiNutrientDetritalFlows (var, reR6c)
!#endif
    end do

end do

!#ifndef isWithBFM
do j = 1, ndetritus
    do i = 1, ndetritus
        if (i /= j) then
            ! calculate flows between detritus groups
            flow2detritus(j) =  flow2detritus(j) &
              + (ep_data(detritus_no(i))%DetPassedProp * ep_data(detritus_no(i))%biomass * ep_detfate(i, j))
        end if
    end do
end do
!#endif

!!!!! initialize the export rate of detritus out of the system

if (FirstTime .eqv. .true.) then
    do j = 1, ndetritus
        det_export_rate(j) = (flow2detritus(j) &
          - es_data(detritus_no(j))%EatenOf + ep_data(detritus_no(j))%detritus_import) / biomass(detritus_no(j))
    end do
end if

FirstTime = .false.

end subroutine calculateDetritalFlows
