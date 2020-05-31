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

subroutine derivs (time, biomass, xdot, biomeq, loss, integrate, lat, lon)

#ifdef isWithBFM
use global_mem
  use mem, ONLY: iiHigherTrophicLevels
#else
  use statevartypesecopath, only: RLEN
#endif

use statevartypesecopath, only: ep_data, flow2detritus, det_export_rate
use statevartypesecosim, only: es_data, es_vul, arena, nvars, imonth, &
                         NutFree, NutTot, NutBiom, NutFreeBase, NutMin, &
                         NutrientForce, detritus_no
#ifdef _ForcePrimaryProd_
use statevartypesecosim, only: PrimaryProdForce
#endif

#ifdef _Ecospace_
use statevartypesecospace, only: QperB, M2
#endif

implicit none

! variables inherited from ecopath model
integer, intent(in)           :: integrate(nvars)
real(RLEN), intent(in)        :: biomass(nvars)   ! array of initial cond.
real(RLEN), intent(in)        :: time
integer, optional, intent(in) :: lat, lon

! in-subroutine variables
integer                 :: var, prey, pred ! loop vars: var, prey, predator
integer                 :: j               ! loop vars
integer                 :: vrows, vcols    ! rows & columns of vul. matrix
real(RLEN)              :: Pmult           ! prim. prod. forcing multip.
real(RLEN)              :: b_derivs(nvars) ! copied b_init(nvars)
real(RLEN)              :: consumption     ! consumption of consumers
real(RLEN)              :: M2_predation    ! predation on groups
#ifdef isWithBFM
  real(RLEN)            :: Predation       ! predation on groups
#endif
real(RLEN), intent(out) :: xdot(nvars)     ! results of state equations
real(RLEN), intent(out) :: biomeq(nvars)   ! rate of change for state vars
real(RLEN), intent(out) :: loss(nvars)     ! sinks for state variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef isWithBFM
call HTLGlobalDynamics()
#endif

vrows = size(es_vul, 1)
vcols = size(es_vul, 2)

! set predator abundance measure 
! for predation rate calculations (setpred)
b_derivs = biomass
do var = 1, nvars
    if (biomass(var) < 1.0e-20) then
        b_derivs(var) = 1.0e-20
    end if  
    if (integrate(var) >= 0) then
        es_data(var)%pred = b_derivs(var)
    end if
end do

! nutrient biomass
NutBiom = 0
do var = 1, nvars
    NutBiom = NutBiom + biomass(var)
end do

if (time == 0. .and. imonth == 0) then
    NutFree = NutTot - NutBiom
else
#ifdef _ForceNutrient_
    if (time == 1.0) then
        NutFree = NutTot - NutBiom
    else
        NutFree = NutTot * NutrientForce(imonth) - NutBiom
    end if 
#endif
#ifndef _ForceNutrient_
    NutFree = NutTot - NutBiom
#endif
end if

! amount of free nutrients (NutFree) must not be less than 
! the background nutrient concentration (NutMin)
if (NutFree < NutMin) then
    NutFree = NutMin
end if

j = 0
do var = 1, nvars
    if (ep_data(var)%org_type /= 0) then
! P/B scaled by Max. Relative Production, i.e. es_data%pob_max.
        if (time == 0. .and. imonth == 0) then
            Pmult = 1.0
            if (ep_data(var)%org_type /= 0) then
                es_data(var)%abs_PoB_max = es_data(var)%rel_PoB_max &
                  * ep_data(var)%PoB
                es_data(var)%PoB_biomass = 0
            end if            
        else
#ifdef _ForcePrimaryProd_
            if (ep_data(var)%org_type == 1 .and. time /= 1.0) then
                j = j + 1
                Pmult = PrimaryProdForce(imonth, j)
            end if
#endif
#ifndef _ForcePrimaryProd_
    Pmult = 1.0
#endif
        end if
        es_data(var)%PoB_base = 2 * NutFree &
          / (NutFree + NutFreeBase(var)) * Pmult &
          * es_data(var)%abs_PoB_max &
         / (1 + biomass(var) * es_data(var)%PoB_biomass)
    end if
end do

! set relative prey switching
 call setRelativeSwitchingParameters (vrows, vcols, biomass)

! vulnerability calculations
 call calculateVulnerableBiomasses (vrows, vcols, biomass)

!!!!! calculate processes; production and consumption
! primary production
do var = 1, nvars
    if (ep_data(var)%org_type == 1) then
        es_data(var)%pp = es_data(var)%PoB_base * biomass(var)
    else
        es_data(var)%pp = 0
    end if
end do

! consumption
do pred = 1, nvars
    if (pred <= vcols) then
        consumption = 0
        do prey = 1, vrows
            if (es_vul(prey, pred) /= -999) then
                consumption = consumption + (arena%a_eff(prey, pred) &
                  * arena%v_biom(prey, pred) * es_data(pred)%pred &
                  / es_data(pred)%hden)
            end if
        end do
        es_data(pred)%qq = consumption + es_data(pred)%QBoutside * biomass(pred)
        es_data(pred)%EatenBy = es_data(pred)%qq
#ifdef _Ecospace_
        if (lat > 0 .and. lon > 0) then
            QperB(lat, lon, pred) = es_data(pred)%EatenBy / biomass(pred)
        end if
#endif
    else
        es_data(pred)%qq = 0
        es_data(pred)%EatenBy = es_data(pred)%qq / biomass(pred)
#ifdef _Ecospace_
        if (lat > 0 .and. lon > 0) then
            QperB(lat, lon, pred) = 0
        end if
#endif
    end if
end do

! unassimilated consumption
do var = 1, nvars
    if (es_data(var)%qq /= 0) then
        es_data(var)%unassimilated = es_data(var)%qq * ep_data(var)%unass_Q
    else
        es_data(var)%unassimilated = 0
    end if
end do

! non-predation mortalities
do var = 1, nvars
    if (ep_data(var)%org_type == 1) then
        es_data(var)%M0 = ((1 - ep_data(var)%EE) * ep_data(var)%PoB) * biomass(var)
    else if (ep_data(var)%org_type == 2) then
        es_data(var)%M0 = ((1 - ep_data(var)%EE) * ep_data(var)%PoB) &
          * (1 - es_data(var)%M0_pred + es_data(var)%M0_pred &
          * es_data(var)%Ftime) * biomass(var) 
    else
        es_data(var)%M0 = 0
    end if
end do

! predation mortalities
do prey = 1, nvars
    M2_predation = 0
#ifdef isWithBFM
    do pred = 1, iiHigherTrophicLevels ! only HTL predators are accounted
#else
    do pred = 1, vcols
#endif
        if (es_vul(prey, pred) /= -999) then
            M2_predation = M2_predation + (arena%a_eff(prey, pred) &
              * arena%v_biom(prey, pred) * es_data(pred)%pred / es_data(pred)%hden)
#ifdef isWithBFM
    Predation = (arena%a_eff(prey, pred) &
        * arena%v_biom(prey, pred) * es_data(pred)%pred / es_data(pred)%hden)
    call calculateCompensatoryFlows (prey, pred, Predation)
#endif
        end if
    end do
    es_data(prey)%M2 = M2_predation
    es_data(prey)%EatenOf = M2_predation
#ifdef _Ecospace_
    if (lat > 0 .and. lon > 0) then
        M2(lat, lon, prey) = es_data(prey)%EatenOf / biomass(prey)
    end if
#endif
end do

! detritus calculations
 call calculateDetritalFlows (biomass)

! sum of sinks for state variables
j = 0
do var = 1, nvars
    if (ep_data(var)%org_type /= 0) then
        loss(var) = es_data(var)%M0 + es_data(var)%M2 + es_data(var)%fishmort &
          * biomass(var)
    else
        j = j + 1
        loss(var) = det_export_rate(j) * biomass(var) + es_data(var)%EatenOf
    end if
end do

! now calculate state equations
j = 0
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        xdot(var) = (ep_data(var)%PoQ * es_data(var)%qq + biomass(var) &
          * es_data(var)%PoB_base) - loss(var)
        if (loss(var) > 0 .and. biomass(var) > 0) then
            biomeq(var) = (ep_data(var)%PoQ * es_data(var)%qq + biomass(var) &
              * es_data(var)%PoB_base) / (loss(var) / biomass(var))
        else
            biomeq(var) = 1.0E-20
        end if
    else if (ep_data(var)%org_type == 1) then
        xdot(var) = es_data(var)%pp - loss(var)
        if (loss(var) > 0 .and. biomass(var) > 0) then
            biomeq(var) = es_data(var)%pp / (loss(var) / biomass(var))
        else
            biomeq(var) = 1.0E-20
        end if
    else
        j = j + 1
        xdot(var) = ep_data(detritus_no(j))%detritus_import + flow2detritus(j) - loss(var)
        if (loss(var) /= 0 .and. biomass(var) > 0 .and. &
          flow2detritus(j) > 0) then
            biomeq(var) = (ep_data(detritus_no(j))%detritus_import + flow2detritus(j)) &
              / (loss(var) / biomass(var))
        else
            biomeq(var) = 1.0E-20
        end if
    end if
end do

end subroutine derivs
