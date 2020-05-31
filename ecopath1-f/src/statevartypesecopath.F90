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

module statevartypesecopath

! This module contains variable definitions and their attributes
! for state variables and statistical objects in the Ecopath model

implicit none

!integer, parameter :: RLEN = selected_real_kind(6,37)
!integer, parameter :: RLEN = selected_real_kind(8,60)
integer, parameter :: RLEN = selected_real_kind(15,307)
!integer, parameter :: RLEN = selected_real_kind(33,4931)
!#ifndef _isDP_
!  integer, parameter :: RLEN = 4
!#endif
!#ifdef _isDP_
!  integer, parameter :: RLEN = 8
!#endif

type ecopath_data

 real(RLEN) :: biomass
 real(RLEN) :: PoB                  ! production/biomass (P/B)
 real(RLEN) :: QoB                  ! consumption/biomass (Q/B)
 real(RLEN) :: EE                   ! ecotrophic efficiency (EE)
 real(RLEN) :: PoQ                  ! production/consumption (P/Q)
 real(RLEN) :: unass_Q              ! unassimilated/consumption
 real(RLEN) :: detritus_import      ! import from outside of the system
 real(RLEN) :: landings             ! landed amount of catch
 real(RLEN) :: discards             ! discarded amount of catch

 integer :: org_type             ! organism type: (1) producer;
                                 ! (2) consumer; (0) detritus

 integer :: isstanza             ! multistanza group or not,
                                 ! (1) multistanza; (0) no stanzas

 integer :: stanza_no            ! multistanza group number
 integer :: age_start            ! age to start stanza (in months)

 integer :: isleading            ! leading stanza or not,
                                 ! (1) leading; (0) substanza

 real(RLEN) :: production           ! B * P/B
 real(RLEN) :: consumption          ! B * Q/B
 real(RLEN) :: respiration          ! (1 - unass_q) * Q - P
 real(RLEN) :: assimilation         ! (1 - unass_q) * Q/B * B
 real(RLEN) :: EatenOf              ! the predation pressure on group
 real(RLEN) :: EatenBy              ! the consumption of group
 real(RLEN) :: DetPassedProp        ! detritus passed on by proportion
 real(RLEN) :: BA                   ! biomass accumulation

end type ecopath_data


type growth_data

integer :: stanza_no
real(RLEN) :: K
real(RLEN) :: RecPow
real(RLEN) :: BaB
real(RLEN) :: WmatWinf

end type growth_data


type multi_stanza

 integer :: ep_groupno(10) ! group numbers of substanzas in multistanza
 real(RLEN) :: biomass(10)
 integer :: age_start(10)  ! Age to start stanza (in months)
 real(RLEN) :: mortality(10)  ! total mortality (Z) of stanzas
 real(RLEN) :: QoB(10)        ! consumption/biomass (Q/B)
 integer :: isleading(10)  ! leading stanza (1) or substanza (0)
 integer :: age_infinity   ! maximum age the leading stanza may attain
 integer :: substanzas     ! number of substanzas in the multistanza
 real(RLEN) :: vbK            ! von Bertalanffy growth coefficient (K)
 real(RLEN) :: vbM            ! von Bertalanffy metabolic parameter
 real(RLEN) :: rec_power      ! recruitment power
 real(RLEN) :: rel_BA         ! relative Biomass Accumulation rate
 real(RLEN) :: Wmat_Winf      ! Ratio of Weight maturity to Weight infinity
 real(RLEN) :: RzeroS         ! Base recruitment to age 0 for split species
 real(RLEN) :: Rhat
 real(RLEN) :: Ahat
 real(RLEN) :: RhatC
 real(RLEN) :: AhatC
 real(RLEN) :: Wage(0:299)    ! relative body weight at age "a" (monthly)
 real(RLEN) :: WWa(0:299)     ! relative consumption (Q) depending on Wage
 real(RLEN) :: survive(0:299) ! survival rate for each age group (monthly)
 real(RLEN) :: splitno(0:299) ! number of survivors at age (monthly)

end type multi_stanza


character ( len = 250 ) :: InputData_fname
character ( len = 250 ) :: DietComp_fname
character ( len = 250 ) :: DetFate_fname
character ( len = 250 ) :: GrowthParam_fname
character ( len = 250 ) :: Results_fname
character ( len = 250 ) :: HDF5_fname
character(len = 250), allocatable :: groupnames(:)
logical                 :: isASCIIinputFile

! flows to detritus (sum of natural mortalities + unass. consumption)
real(RLEN), allocatable :: m_flows2det(:,:)
! matrix of predations for each group
real(RLEN), allocatable :: m_consumed(:)

! Ecopath data matrix
type(ecopath_data), allocatable, target :: ep_data(:)

! Multistanza growth parameter matrix
type(growth_data), allocatable, target  :: ep_growth(:)

! data matrix of multistanza parameters and calculations
type(multi_stanza), allocatable, target :: ms_data(:)
real(RLEN), allocatable, target :: ep_diet(:, :) ! Ecopath diet matrix
real(RLEN), allocatable, target :: ep_detfate(:, :) ! Ecopath detritus fate matrix

integer              :: nvars, nstanzas

end module statevartypesecopath
