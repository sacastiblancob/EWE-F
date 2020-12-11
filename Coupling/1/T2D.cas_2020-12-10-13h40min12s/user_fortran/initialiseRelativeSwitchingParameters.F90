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

         subroutine initialiseRelativeSwitchingParameters &
  (nvars, vrows, vcols, ep_data, es_data, es_vul, arena)

!use statevartypesecopath, only: ep_data
use statevartypesecopath, only: RLEN, ecopath_data
!use statevartypesecosim, only: es_data, es_vul, arena
use statevartypesecosim, only: ecosim_data, arena_data

implicit none

! INTENT VARIABLES
integer, intent(in) :: nvars, vrows, vcols
TYPE(ecopath_data), INTENT(IN)  :: ep_data(nvars)
TYPE(ecosim_data), INTENT(INOUT) :: es_data(nvars)
REAL(RLEN), INTENT(IN) :: es_vul(vrows, vcols)
TYPE(arena_data), INTENT(INOUT) :: arena

! IN SUBROUTINE VARIABLES
integer             :: prey, pred

!!!!! initial switching parameters (InitRelaSwitch)
do pred = 1, vcols
    es_data(pred)%pred_den = 0
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            es_data(pred)%pred_den = es_data(pred)%pred_den &
              + arena%a(prey, pred) * ep_data(prey)%biomass &
              ** es_data(pred)%switch_power
        end if
    end do
end do

do pred = 1, vcols
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
!#ifdef isWithBFM
!            arena%base_time_switch(prey, pred) = arena%a(prey, pred) &
!              * ep_data(prey)%biomass ** es_data(pred)%switch_power &
!              / (es_data(pred)%pred_den + 1.0e-20)
!#else
!            arena%base_time_switch(prey, pred) = real(arena%a(prey, pred) &
!              * ep_data(prey)%biomass ** es_data(pred)%switch_power &
!              / (es_data(pred)%pred_den + 1.0e-20), 4)
            arena%base_time_switch(prey, pred) = real(arena%a(prey, pred) &
              * ep_data(prey)%biomass ** es_data(pred)%switch_power &
              / (es_data(pred)%pred_den + 1.0e-20), RLEN)
!#endif
        end if
    end do
end do

end subroutine initialiseRelativeSwitchingParameters
