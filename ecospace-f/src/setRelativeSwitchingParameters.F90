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

subroutine setRelativeSwitchingParameters (vrows, vcols, biomass)

#ifdef isWithBFM
  use global_mem
#else
  use statevartypesecopath, only: RLEN
#endif

use statevartypesecosim, only: es_data, es_vul, arena, nvars

implicit none

integer, intent(in) :: vrows, vcols
real(RLEN), intent(in) :: biomass(nvars)
integer             :: prey, pred

do pred = 1, vcols
    es_data(pred)%pred_den = 0
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            es_data(pred)%pred_den = es_data(pred)%pred_den + arena%a(prey, pred) &
              * biomass(prey) ** es_data(pred)%switch_power
            arena%rela_switch(prey, pred) = 1
        end if
    end do
end do

do pred = 1, vcols
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            if (es_data(pred)%switch_power > 0) then
#ifdef isWithBFM
                arena%rela_switch(prey, pred) = arena%a(prey, pred) &
                  * biomass(prey) ** es_data(pred)%switch_power &
                  / (es_data(pred)%pred_den + 1.0e-20) &
                  / arena%base_time_switch(prey, pred)
#else
                arena%rela_switch(prey, pred) = real(arena%a(prey, pred) &
                  * biomass(prey) ** es_data(pred)%switch_power &
                  / (es_data(pred)%pred_den + 1.0e-20) &
                  / arena%base_time_switch(prey, pred), 4)
#endif 
            end if
        end if
    end do
end do


end subroutine setRelativeSwitchingParameters
