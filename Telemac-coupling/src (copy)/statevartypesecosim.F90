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

module statevartypesecosim

! This module contains variable definitions and their attributes
! for state variables and foragin arena parameters in the Ecosim model.

!#ifdef isWithBFM
!  use global_mem
!#endif

!#ifndef isWithBFM
  use statevartypesecopath, only: RLEN
!#endif

implicit none

type ecosim_data

! parts inherited from scenario parameters input file
 real(RLEN) :: rel_PoB_max     ! maximum relative P/B
 real(RLEN) :: Ftime           ! feeding time
 real(RLEN) :: Ftime_max       ! maximum relative feeding time
 real(RLEN) :: Ftime_adjust    ! feeding time adjustment rate

 real(RLEN) :: M0_pred         ! fraction of mortality sensitive
                            ! to changes in feeding time

 real(RLEN) :: risk_time       ! predator effect on feeding time
 real(RLEN) :: Q_maxoQ_0       ! density dependent catchability
 real(RLEN) :: QB_maxoQB_0     ! maximum relative consumption
 real(RLEN) :: switch_power    ! prey switching power parameter
 integer    :: isAdvected      ! if the group is impacted by advective flows

! parts that are calculated within the Ecosim model
 real(RLEN) :: EatenOf         ! the predation pressure on group in time
 real(RLEN) :: EatenBy         ! the consumption of group in time

 real(RLEN) :: CB_base         ! base consumption biomass ratio
                            ! calculated from initial conditions

 real(RLEN) :: CB_last         ! last consumption biomass ratio
 real(RLEN) :: QBoutside       ! QoB value for groups with import in diets
 real(RLEN) :: Q_main
 real(RLEN) :: Q_risk
 real(RLEN) :: Q_opt
 real(RLEN) :: risk_rate
 real(8)    :: pred_den        ! Denominator for calculating
                            ! prey switching parameters

 real(RLEN) :: hden            ! Actually QBmaxQBo/(QbmaxQBo-1);
                            ! but at time=0 => QBmaxQBo/(QbmaxQBo+1)
 real(RLEN) :: hdent

 real(RLEN) :: PoB_base        ! actual P/B for primary production;
                            ! scaled by nutrient concentration

 real(RLEN) :: PoB_biomass     ! P/B scaled by producer biomass
 real(RLEN) :: abs_PoB_max     ! maximum absolute P/B
 real(RLEN) :: htime           ! Handling time for predators

! sources and sinks summed
 real(RLEN) :: pp              ! primary production for producers
 real(RLEN) :: M2              ! predation mortality of state variables
 real(RLEN) :: M0              ! non-predation natural mortality
 real(RLEN) :: qq              ! consumption of consumers
 real(RLEN) :: unassimilated   ! unassimilated ration of consumed food
 real(RLEN) :: pred

! parameters for calculating fishing mortality
 real(RLEN) :: Q_mult          ! multiplier: density-dependent catchability
 real(RLEN) :: fishmort        ! fishing mortality

end type ecosim_data



type ecosim_multi_stanza

 real(RLEN) :: NageS(0:299)
 real(RLEN) :: WageS(0:299)
 real(RLEN) :: EggsSplit(0:299)
 real(RLEN) :: SplitAlpha(0:299)
 real(RLEN) :: BaseEggsStanza
 real(RLEN) :: EggsStanza
 real(RLEN) :: RscaleSplit

end type ecosim_multi_stanza



! these are attributes related to specific foraging arenas
type arena_data

! Search rate of predator j for prey i
 real(RLEN), allocatable :: a(:, :)

! Consumption of predator j of prey i in arena n
 real(RLEN), allocatable :: Q_link(:, :)

! Sum of consumptions of predator j in arena n
 real(RLEN), allocatable :: Q_arena(:, :)

! Vulnerability parameter between predator j and prey i in arena n
 real(RLEN), allocatable :: vul_arena(:, :)

! Vulnerable biomass of prey i in arena n to predator j
 real(RLEN), allocatable :: vul_biom(:, :)

! Vulnerability rate calculated from vulnerability multiplier matrix
 real(RLEN), allocatable :: vulrate(:, :)

! Search rate of predator j for prey i in arena n
 real(RLEN), allocatable :: a_link(:, :)

! Basic prey switching time parameter for predator j
 real(RLEN), allocatable :: base_time_switch(:, :)

! Relative prey switching parameter for predator j
 real(RLEN), allocatable :: rela_switch(:, :)

! Effective search rate of predator j
 real(RLEN), allocatable :: a_eff(:, :)

! Effective vulnerability rate for predator j
 real(RLEN), allocatable :: v_eff(:, :)

! Actual vulnerable biomass of prey i to predator j
 real(RLEN), allocatable :: v_biom(:, :)

! Denominator for calculating vulnerable biomass
 real(RLEN), allocatable :: v_denom(:, :)

end type arena_data



type forcing_data

! time series of forced biomass values
real(RLEN), allocatable :: biomass(:, :)

! time series of forced catch values
real(RLEN), allocatable :: catches(:, :)

! time series of forced fishing mortality (F) values
real(RLEN), allocatable :: fishforce(:, :)

! type of forcing data
integer, allocatable :: forcetype(:)

end type forcing_data



! Ecosim data matrix (see statevartypesecosim.f90)
type(ecosim_data), allocatable         :: es_data(:)

type(ecosim_multi_stanza), allocatable :: es_ms_data(:)

! Ecosim vulnerability matrix, read from file
real(RLEN), allocatable                   :: es_vul(:, :)

! parameters for foraging arenas (see statevartypesecosim.f90)
type(arena_data)                       :: arena

! time series of forcing data; fishing mortality(F)
type(forcing_data)                     :: force

logical                                :: FirstTime, UpdateStanzas

! nutrients
real(RLEN) :: NutTot                      ! total amount of nutrients
real(RLEN) :: NutBiom                     ! biomass of nutrients
real(RLEN) :: NutFree                     ! amount of free nutrients

! background concentration for nutrients
real(RLEN) :: NutMin

! base concentration of free nutrients
real(RLEN), allocatable :: NutFreeBase(:)

real(RLEN), allocatable :: NutrientForce(:)

real(RLEN), allocatable :: PrimaryProdForce(:, :)

! nutrient concentration inputs
real(RLEN) :: NutBaseFreeProp  ! base proportion of free nutrients
real(RLEN) :: NutPBmax         ! max P/B due to nutrient concentration
integer :: StepsPerMonth    ! number of time steps per month
real(RLEN) :: relaxeco            ! relaxation parameter
!tf original is integer
real(RLEN) :: tf
!integer :: tf               ! number of years to simulate
integer :: imonth           ! integer month

character(len = 250) :: GroupInfo_fname
character(len = 250) :: Vulnerability_fname
character(len = 250) :: Forcing_fname
character(len = 250) :: NutrientForcing_fname
character(len = 250) :: PrimaryProdForcing_fname
character(len = 250), allocatable :: groupnames(:)

integer, allocatable :: detritus_no(:)
integer :: ndetritus
integer, allocatable :: producer_no(:)
integer :: nproducer

real(RLEN), allocatable :: BB(:)
real(RLEN), allocatable :: BBAvg(:), LossAvg(:)
real(RLEN), allocatable :: EatenByAvg(:)
real(RLEN), allocatable :: EatenOfAvg(:)
real(RLEN), allocatable :: PredAvg(:)

integer              :: nvars, nstanzas

! logicals for know if compute Primary Production Forcing and Nutrient Forcing
  logical    :: boolFN = .false.
  logical    :: boolFPP = .false.

!!!!ALL NEXT VARIABLES COMES FROM ECOSIM.F90 ORIGINAL FILE, DEFINED DOWN IMPLICIT NONE

! loops vars
  integer              :: iec, j, m, n, var  ! loop vars; i prey & j predator

! initial time, number of time steps, size of one step, time step and time
  integer              :: t0               ! simulation start time
  integer              :: noftsteps, step  ! number of time steps
  real(RLEN)           :: tstep            ! time step
  real(RLEN)           :: time
  real(RLEN)           :: SecondsPerMonth  !! Second Per Month for Telemac Coupling

!  ! (1) integrate, (0) do not integrate
  integer, allocatable :: integrate(:)

  ! intital biomass values
  real(RLEN), allocatable :: B(:), b_pred(:)  ! initial biomass values

  ! array of changes in biomass values calculated by derivs()
  real(RLEN), allocatable :: xdot(:)

  ! array of changes in non-integrated biomass values calculated by derivs
  real(RLEN), allocatable :: biomeq(:)

  ! array of sinks for state variables calculated by derivs
  real(RLEN), allocatable :: loss(:)

  ! array of integrated biomass values calculated by rk4 in each time step
  real(RLEN), allocatable :: b_out(:)

  ! rate of change in state variables within the time step
  real(RLEN), allocatable :: rrate(:)

!#ifdef isWithBFM
!  real(RLEN)              ::  multistanza_update
!  real(RLEN), allocatable ::  ruHTLc(:,:,:),ruHTLn(:,:,:),ruHTLp(:,:,:)
!  real(RLEN), allocatable ::  ruHTLl(:,:,:),ruHTLi(:,:,:),ruHTLs(:,:,:)
!  real(RLEN), allocatable ::  tfluxc(:,:), tfluxn(:,:), tfluxp(:,:)
!  real(RLEN), allocatable ::  p_qncHTL(:,:),  p_qpcHTL(:,:)
!  real(RLEN), allocatable ::  flow2detritusR6c(:,:), flow2detritusR6n(:,:), flow2detritusR6p(:,:)
!#endif

end module statevartypesecosim
