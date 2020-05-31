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

subroutine calculateLotkaVolterraEffectiveSearchRates (vrows, vcols)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

use statevartypesecopath, only: ep_data, ep_diet
use statevartypesecosim, only: es_data, es_vul, arena, nvars

implicit none

integer, intent(in) :: vrows, vcols
integer             :: var, prey, pred
real(RLEN)             :: Dzero, Denv, Consumption

do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
!#ifdef isWithBFM
!        es_data(var)%htime = es_data(var)%pred &
!          / (es_data(var)%QB_maxoQB_0 &
!          * ep_data(var)%biomass * ep_data(var)%QoB)
!#else
        es_data(var)%htime = real(dble(es_data(var)%pred) &
          / (dble(es_data(var)%QB_maxoQB_0) &
          * dble(ep_data(var)%biomass) * dble(ep_data(var)%QoB)), 4)
!#endif
    else
        es_data(var)%htime = 0
    end if
end do

do pred = 1, vcols
    do prey = 1, vrows

        Dzero = es_data(pred)%QB_maxoQB_0 / (es_data(pred)%QB_maxoQB_0 - 1)

        if (es_vul(prey, pred) /= -999) then
            Consumption = (ep_data(pred)%biomass * ep_data(pred)%QoB &
              * ep_diet(prey, pred))
!#ifdef isWithBFM
!            arena%vulrate(prey, pred) = es_vul(prey, pred) &
!              * Consumption / ( ep_data(prey)%biomass+ 0.0001 )
!#else
            arena%vulrate(prey, pred) = real(dble(es_vul(prey, pred)) &
                 * dble(Consumption) / dble(ep_data(prey)%biomass), 4)
!#endif

            ! Denominator for calculating search rate of predator
            Denv = ep_data(prey)%biomass &
              * es_data(pred)%pred * arena%vulrate(prey, pred) &
              - Consumption * es_data(pred)%pred

            if (Denv < 1.0e-20) then
                Denv = 1.0e-20
            end if

            arena%a(prey, pred) = Dzero * 2 &
              * Consumption * arena%vulrate(prey, pred) &
              / Denv
        else
            ! arena%vulrate(prey, pred) = 0
            arena%a(prey, pred) = 0
        end if
    end do
end do

end subroutine calculateLotkaVolterraEffectiveSearchRates
