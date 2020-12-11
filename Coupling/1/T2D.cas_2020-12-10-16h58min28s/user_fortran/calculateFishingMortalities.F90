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

            subroutine calculateFishingMortalities &
  (nvars, imonth, biomass, force, ep_data, es_data)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

!use statevartypesecopath, only: ep_data
use statevartypesecopath, only: ecopath_data
!use statevartypesecosim, only: es_data, force, nvars, imonth
use statevartypesecosim, only: forcing_data, ecosim_data

implicit none

! INTENT VARIABLES
INTEGER, INTENT(IN) :: nvars, imonth
REAL(RLEN), intent(in) :: biomass(nvars)
TYPE(forcing_data), INTENT(IN) :: force
TYPE(ecopath_data), INTENT(IN) :: ep_data(nvars)
TYPE(ecosim_data), INTENT(INOUT) :: es_data(nvars)

! IN SUBROUTINE VARIABLES
integer             :: var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do var = 1, nvars
    if (es_data(var)%Q_maxoQ_0 /= -999) then
        es_data(var)%Q_mult = es_data(var)%Q_maxoQ_0 &
          / (1 + (es_data(var)%Q_maxoQ_0 - 1) &
          * biomass(var) / ep_data(var)%biomass)
    else
        es_data(var)%Q_mult = 0
    end if
    if (force%forcetype(var) == 5 .and. force%fishforce(imonth, var) /= -999) then
        es_data(var)%fishmort = ((ep_data(var)%landings + ep_data(var)%discards) &
          / ep_data(var)%biomass) * es_data(var)%Q_mult * force%fishforce(imonth, var)
    else if (force%forcetype(var) == 4 .and. force%fishforce(imonth, var) /= -999) then
        es_data(var)%fishmort = force%fishforce(imonth, var) &
          * es_data(var)%Q_mult
    else if  (force%fishforce(imonth, var) == -999) then
        es_data(var)%fishmort = ((ep_data(var)%landings + ep_data(var)%discards) &
          / ep_data(var)%biomass) * es_data(var)%Q_mult
    end if
end do

end subroutine calculateFishingMortalities
