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

              subroutine calculateNutrientConcentrations &
          (nvars, ep_data, NutFreeBase, NutBaseFreeProp, NutMin, &
           NutTot, NutBiom, NutFree)

!use statevartypesecopath, only: ep_data
use statevartypesecopath, only: RLEN, ecopath_data
!use statevartypesecosim, only: nvars, NutBaseFreeProp, &
!                         NutFreeBase, NutMin, NutTot, NutBiom, NutFree
!use statevartypesecosim, only: NutBaseFreeProp, &
!                         NutFreeBase, NutMin, NutTot, NutBiom, NutFree

implicit none

! INTENT VARIABLES
INTEGER, INTENT(IN)             :: nvars
TYPE(ecopath_data), INTENT(IN)  :: ep_data(nvars)
REAL(RLEN), INTENT(INOUT)       :: NutFreeBase(nvars)
REAL(RLEN), INTENT(INOUT)       :: NutBiom, NutTot, NutFree, &
                                   NutBaseFreeProp, NutMin

! IN SUBROUTINE VARIABLES
integer             :: var

! calculate nutrient biomass
NutBiom = 0
do var = 1, nvars
    NutBiom = NutBiom + ep_data(var)%biomass
!    write(*,*) NutBiom
end do

!write(*,*)
!write(*,*) NutBaseFreeProp

NutTot  = NutBiom / (1 - NutBaseFreeProp)
NutFree = NutTot - NutBiom

!write(*,*) NutTot
!write(*,*) NutFree

! base concentration of free nutrients
!allocate(NutFreeBase(nvars))
do var = 1, nvars
    if (ep_data(var)%org_type /= 0) then
        ! NutFreeBase(var) = (es_data(var)%rel_PoB_max - 1) * NutFree
        NutFreeBase(var) = NutFree
    end if
end do

!#ifdef isWithBFM
!  NutMin = 0.00101D0 * NutFree
!#else
  NutMin = real(0.00101D0 * NutFree, RLEN)
!#endif

end subroutine calculateNutrientConcentrations
