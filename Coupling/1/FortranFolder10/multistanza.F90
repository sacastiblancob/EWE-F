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

                          subroutine multistanza &
  (nvars, nstanzas, BtoUpdate, ms_data, es_data, es_ms_data, BBAvg, LossAvg)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

!use statevartypesecopath, only: ms_data
  use statevartypesecopath, only: multi_stanza
!use statevartypesecosim, only: es_data, es_ms_data, nvars, nstanzas, &
!                         BBAvg, LossAvg
  use statevartypesecosim, only: ecosim_data, ecosim_multi_stanza

implicit none

! variables inherited
INTEGER, INTENT(IN)                      :: nvars, nstanzas
real(RLEN), intent(inout)                :: BtoUpdate(nvars)
TYPE(multi_stanza), INTENT(IN)           :: ms_data(nstanzas)
TYPE(ecosim_data), INTENT(IN)            :: es_data(nvars)
TYPE(ecosim_multi_stanza), INTENT(INOUT) :: es_ms_data(nstanzas)
REAL(RLEN), INTENT(IN), DIMENSION(nvars) :: BBAvg, LossAvg


! in-subroutine variables
integer :: substanza         ! substanza number
integer :: AgeMax            ! maximum age of substanza
integer :: AgeMin            ! minimum age of substanza
integer :: age               ! age (in months) of substanzas
integer :: stanza            ! stanza number
integer :: age_last          ! maximum age of substanza
real(RLEN) :: Be, Gf, Su        ! what are these??
real(RLEN) :: Nt                ! numbers of animals in substanza


do stanza = 1, nstanzas

    Be = 0
    do substanza = 1, ms_data(stanza)%substanzas

!#ifdef isWithBFM
!        Su = exp(-LossAvg(ms_data(stanza)%ep_groupno(substanza)) / 12.0D0&
!          / BBAvg(ms_data(stanza)%ep_groupno(substanza)))
!#else
        Su = real(exp(-LossAvg(ms_data(stanza)%ep_groupno(substanza)) / 12.0D0&
          / BBAvg(ms_data(stanza)%ep_groupno(substanza))), 4)
!#endif

        Gf = es_data(ms_data(stanza)%ep_groupno(substanza))%EatenBy &
          / es_data(ms_data(stanza)%ep_groupno(substanza))%pred

        if (substanza < ms_data(stanza)%substanzas) then
            age_last = ms_data(stanza)%age_start(substanza+1) - 1
        else
            age_last = ms_data(stanza)%age_infinity
        end if

        do age = ms_data(stanza)%age_start(substanza), age_last
            es_ms_data(stanza)%NageS(age) &
              = es_ms_data(stanza)%NageS(age) * Su
            es_ms_data(stanza)%WageS(age) = ms_data(stanza)%vbM &
              * es_ms_data(stanza)%WageS(age) &
              + Gf * es_ms_data(stanza)%SplitAlpha(age)

            if (es_ms_data(stanza)%WageS(age) &
              > ms_data(stanza)%Wmat_Winf) then
                Be = Be + es_ms_data(stanza)%NageS(age) * &
                  (es_ms_data(stanza)%WageS(age) &
                  - ms_data(stanza)%Wmat_Winf)
            end if
        end do
    end do

    es_ms_data(stanza)%WageS(ms_data(stanza)%age_infinity) &
      = (Su * ms_data(stanza)%Ahat &
      + (1 - Su) &
      * es_ms_data(stanza)%WageS(ms_data(stanza)%age_infinity - 1)) &
      / (1 - ms_data(stanza)%Rhat * Su)

    es_ms_data(stanza)%EggsStanza = Be

    do substanza = ms_data(stanza)%substanzas, 1, -1

        if (substanza == ms_data(stanza)%substanzas) then
            AgeMax = ms_data(stanza)%age_infinity
        else
            AgeMax = ms_data(stanza)%age_start(substanza + 1) - 1
        end if

        if (substanza > 1) then
            AgeMin = ms_data(stanza)%age_start(substanza)
        else
            AgeMin = 1
        end if

        if (substanza == ms_data(stanza)%substanzas) then
            Nt = es_ms_data(stanza)%NageS(AgeMax) &
              + es_ms_data(stanza)%NageS(AgeMax - 1)

            if (Nt == 0) then
                Nt = 1.0e-30
            end if

            es_ms_data(stanza)%NageS(AgeMax) = Nt
            AgeMax = AgeMax - 1
        end if

        do age = AgeMax, AgeMin, -1
            es_ms_data(stanza)%NageS(age) &
              = es_ms_data(stanza)%NageS(age - 1)
            es_ms_data(stanza)%WageS(age) &
              = es_ms_data(stanza)%WageS(age - 1)
        end do

    end do

    if (es_ms_data(stanza)%BaseEggsStanza > 0) then
        es_ms_data(stanza)%NageS(ms_data(stanza)%age_start(1)) &
          = es_ms_data(stanza)%RscaleSplit &
          * 1. * ms_data(stanza)%RzeroS * 1.
    end if

!#ifdef isWithBFM
!    es_ms_data(stanza)%NageS(ms_data(stanza)%age_start(1)) = &
!      es_ms_data(stanza)%NageS(ms_data(stanza)%age_start(1)) &
!      * es_ms_data(stanza)%EggsStanza &
!      / es_ms_data(stanza)%BaseEggsStanza &
!      ** ms_data(stanza)%rec_power
!#else
    es_ms_data(stanza)%NageS(ms_data(stanza)%age_start(1)) = &
      real(dble(es_ms_data(stanza)%NageS(ms_data(stanza)%age_start(1))) &
      * (dble(es_ms_data(stanza)%EggsStanza) &
      / dble(es_ms_data(stanza)%BaseEggsStanza)) &
      ** dble(ms_data(stanza)%rec_power), 4)
!#endif

    es_ms_data(stanza)%WageS(ms_data(stanza)%age_start(1)) = 0
end do

 call setSplitPred(nvars, nstanzas, BtoUpdate, ms_data, es_data, es_ms_data)

end subroutine multistanza
