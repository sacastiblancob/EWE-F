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

subroutine initialiseSplitGroups (B)

#ifdef isWithBFM
  use global_mem
#else
  use statevartypesecopath, only: RLEN
#endif

! This subroutine initialize the parameters of the multistanza groups.

use statevartypesecopath, only: ep_data, ms_data
use statevartypesecosim, only: es_data, es_ms_data, nvars, nstanzas

implicit none

! variables inherited from Ecopath-FORTRAN
real(RLEN), intent(inout) :: B(nvars)

! in-subrotuine variables
real(RLEN) :: Be         ! what's this?
integer :: age        ! age of each substanza in multistanza group
integer :: Agem       ! max. age of each substanza in multistanza group
integer :: stanza     ! stanza number of multistanza group
integer :: substanza  ! substanza number of each substanza

do stanza = 1, nstanzas

    Be = 0
    do age = 0, ms_data(stanza)%age_infinity
    
        es_ms_data(stanza)%NageS(age) = ms_data(stanza)%SplitNo(age)
        es_ms_data(stanza)%WageS(age) = ms_data(stanza)%Wage(age)
    
        if (es_ms_data(stanza)%WageS(age) &
          > ms_data(stanza)%Wmat_Winf) then

            es_ms_data(stanza)%EggsSplit(age) = &
              es_ms_data(stanza)%WageS(age) - ms_data(stanza)%Wmat_Winf

            Be = Be + es_ms_data(stanza)%NageS(age) &
              * es_ms_data(stanza)%EggsSplit(age)
        end if
    end do

    es_ms_data(stanza)%BaseEggsStanza = Be
    es_ms_data(stanza)%EggsStanza     = Be
    es_ms_data(stanza)%RscaleSplit    = 1
end do

! calculate the number of biomasses in split groups
 call setSplitPred (B)

do stanza = 1, nstanzas

    do substanza = 1, ms_data(stanza)%substanzas

        if (substanza < ms_data(stanza)%substanzas) then
            Agem = ms_data(stanza)%age_start(substanza + 1) - 1
        else
            Agem = ms_data(stanza)%age_infinity
        end if
    
        if (substanza == ms_data(stanza)%substanzas) then
            Agem = ms_data(stanza)%age_infinity - 1
        end if

        do age = ms_data(stanza)%age_start(1), Agem
            es_ms_data(stanza)%SplitAlpha(age) = &
              (ms_data(stanza)%Wage(age + 1) &
              - ms_data(stanza)%vbM * ms_data(stanza)%Wage(age)) &
              * es_data(ms_data(stanza)%ep_groupno(substanza))%pred &
              / ep_data(ms_data(stanza)%ep_groupno(substanza))%EatenBy
        end do
    end do

    es_ms_data(stanza)%SplitAlpha(ms_data(stanza)%age_infinity ) = &
      es_ms_data(stanza)%SplitAlpha(ms_data(stanza)%age_infinity - 1)

end do

end subroutine initialiseSplitGroups
