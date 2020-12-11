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

                    subroutine updateForagingTimes &
  (nvars, integrate, ep_data, es_data, EatenByAvg, EatenOfAvg, BBAvg, PredAvg)

!use statevartypesecopath, only: ep_data, RLEN
use statevartypesecopath, only: ecopath_data, RLEN
!use statevartypesecosim, only: es_data, nvars, EatenByAvg, &
!                         EatenOfAvg, BBAvg, PredAvg
use statevartypesecosim, only: ecosim_data

implicit none

! INTENT VARIABLES
INTEGER, INTENT(IN) :: nvars
integer, intent(in) :: integrate(nvars)
TYPE(ecopath_data), INTENT(IN) :: ep_data(nvars)
TYPE(ecosim_data), INTENT(INOUT) :: es_data(nvars)
REAL(RLEN), INTENT(IN), DIMENSION(nvars) :: EatenByAvg, EatenOfAvg, BBAvg, PredAvg

! IN SUBROUTINE VARIABLES
integer :: var

do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        if (integrate(var) == var .or. integrate(var) < 0) then
            es_data(var)%CB_last = EatenByAvg(var) / PredAvg(var)
        end if

        if (ep_data(var)%PoB /= -999) then
!#ifdef isWithBFM
!            es_data(var)%risk_rate = EatenOfAvg(var) / BBAvg(var) &
!              + ((1 - ep_data(var)%EE) * ep_data(var)%PoB) + 0.0000000001D0
!#else
!!            es_data(var)%risk_rate = real(EatenOfAvg(var) / BBAvg(var) &
!!             + ((1 - ep_data(var)%EE) * ep_data(var)%PoB) + 0.0000000001D0, 4)
            es_data(var)%risk_rate = real(EatenOfAvg(var) / BBAvg(var) &
              + ((1 - ep_data(var)%EE) * ep_data(var)%PoB) + 0.0000000001D0, RLEN)
!#endif
        else
!#ifdef isWithBFM
!            es_data(var)%risk_rate = 0.0000000001D0
!#else
!!            es_data(var)%risk_rate = real(0.0000000001D0, 4)
             es_data(var)%risk_rate = real(0.0000000001D0, RLEN)
!#endif
        end if

        es_data(var)%Q_opt = es_data(var)%Q_main + es_data(var)%Q_risk &
          / es_data(var)%risk_rate

        if (es_data(var)%CB_last > 0 .and. (integrate(var) == var &
             .or. integrate(var) < 0)) then

!#ifdef isWithBFM
!            es_data(var)%Ftime = 0.1D0 + 0.9D0 * es_data(var)%Ftime &
!              * (1 - es_data(var)%Ftime_adjust &
!              + es_data(var)%Ftime_adjust * es_data(var)%Q_opt &
!              / es_data(var)%CB_last)
!#else
            es_data(var)%Ftime = real(0.1D0 + 0.9D0 * es_data(var)%Ftime &
              * (1 - es_data(var)%Ftime_adjust &
              + es_data(var)%Ftime_adjust * es_data(var)%Q_opt &
              / es_data(var)%CB_last), 4)
!#endif

            if (es_data(var)%Ftime > es_data(var)%Ftime_max) then
                es_data(var)%Ftime = es_data(var)%Ftime_max
            end if
        end if
    end if
end do

end subroutine updateForagingTimes
