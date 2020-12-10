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

subroutine calculateVulnerableBiomasses (vrows, vcols, biomass)

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
real(RLEN)             :: dwe


do pred = 1, vcols
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            es_data(pred)%hdent = 0
            arena%v_denom(prey, pred) = 0
        end if
    end do
end do

dwe = 0.5

do pred = 1, vcols
    do prey = 1, vrows
        arena%v_denom(prey, pred) = 0
        if (es_vul(prey, pred) /= -999) then
            arena%a_eff(prey, pred) = arena%a_link(prey, pred) * es_data(pred)%Ftime &
              * arena%rela_switch(prey, pred)

            arena%v_eff(prey, pred) = arena%vul_arena(prey, pred) * es_data(prey)%Ftime

            arena%v_denom(prey, pred) = arena%v_denom(prey, pred) + arena%a_eff(prey, pred) &
              * es_data(pred)%pred / es_data(pred)%hden
             ! WRITE(*,*) 'arena_a_eff', arena%a_eff(prey, pred)
             ! WRITE(*,*) 'arena_v_eff', arena%v_eff(prey, pred)
             ! WRITE(*,*) 'arena_v_denom', arena%v_denom(prey, pred)
             ! WRITE(*,*) 'es_hden ', es_data(pred)%hden
        end if
    end do
end do

do pred = 1, vcols
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            arena%v_biom(prey, pred) = arena%v_eff(prey, pred) * biomass(prey) &
              / (arena%vul_arena(prey, pred) + arena%v_eff(prey, pred) &
              + arena%v_denom(prey, pred))
             ! WRITE(*,*) 'arena_v_biom', arena%v_biom(prey, pred)
        end if
    end do
end do

do pred = 1, vcols
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            es_data(pred)%hdent = es_data(pred)%hdent + arena%a_eff(prey, pred) &
              * arena%v_biom(prey, pred)
        end if
    end do
    !WRITE(*,*) 'ES_HDENT ', es_data(pred)%hdent
end do

do pred = 1, vcols
#ifdef isWithBFM
    es_data(pred)%hden = (1 - dwe) * (1 + es_data(pred)%htime &
      * es_data(pred)%hdent) + dwe * es_data(pred)%hden
#else
    es_data(pred)%hden = (1 - dwe) * real((1 + dble(es_data(pred)%htime) &
      * dble(es_data(pred)%hdent)),4) + dwe * es_data(pred)%hden
      !WRITE(*,*) 'ES_HDEN ', es_data(pred)%hden
      !WRITE(*,*) 'ES_Htime ', es_data(pred)%htime
      !WRITE(*,*) 'ES_Hdent ', es_data(pred)%hdent
#endif
end do

do pred = 1, vcols
    do prey = 1, vrows
        arena%v_denom(prey, pred) = 0
        if (es_vul(prey, pred) /= -999) then
            arena%v_denom(prey, pred) = arena%v_denom(prey, pred) + arena%a_eff(prey, pred) &
              * es_data(pred)%pred / es_data(pred)%hden
            !WRITE(*,*) 'arena_v_denom ', arena%v_denom(prey, pred)
        end if
    end do
end do

do pred = 1, vcols
    do prey = 1, vrows
        if (es_vul(prey, pred) /= -999) then
            arena%v_biom(prey, pred) = arena%v_eff(prey, pred) * biomass(prey) &
              / (arena%vul_arena(prey, pred) + arena%v_eff(prey, pred) &
              + arena%v_denom(prey, pred))
            !WRITE(*,*) 'arena_v_biom ', arena%v_biom(prey, pred)
        end if
    end do
end do


end subroutine calculateVulnerableBiomasses
