!                    ***************************
                     MODULE DECLARATIONS_ECOPATH
!                    ***************************

  use statevartypesecopath, only: RLEN

  use statevartypesecospace

  ! This is the main ecosim model. All subroutines are called from
  ! this model and state equations are calculated here.

  use statevartypesecopath, only: HDF5_fname, ep_data, ms_data, &
       det_export_rate, flow2detritus
  use statevartypesecosim, only: GroupInfo_fname, &
       Vulnerability_fname, Forcing_fname, &
       NutrientForcing_fname, PrimaryProdForcing_fname, &
       tf, StepsPerMonth, NutBaseFreeProp, NutPBmax, &
       relax, es_data, nvars, nstanzas, ndetritus, &
       detritus_no, BBAvg, &
       EatenOfAvg, EatenByAvg, PredAvg, LossAvg, &
       es_ms_data, imonth, FirstTime, UpdateStanzas, BB, groupnames
  use readHDF5Database

  implicit none

  character(len = 2000):: FMT0, FMT1, FMT2
  integer              :: vrows, vcols     ! rows & columns of vul. matrix
  integer              :: i, j, m, n, var  ! loop vars; i prey & j predator

  integer              :: lat, lon         ! coordinate variables

  integer              :: t0               ! simulation start time
  integer              :: noftsteps, step  ! number of time steps
  real(RLEN)           :: tstep            ! time step
  real(RLEN)           :: time

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

  ! (1) integrate, (0) do not integrate
  integer, allocatable :: integrate(:)

!ifdef _Ecospace_
  ! matrix of biomass results absolute
  real(RLEN), allocatable :: mat_out(:, :, :, :)

  ! matrix of biomass results relative
  real(RLEN), allocatable :: rel_out(:, :, :, :)

  ! matrix of monthly biomass and catch results absolute
  real(RLEN), allocatable :: mat_out_monthly (:, :, :, :)
  real(RLEN), allocatable :: catch_out_monthly (:, :, :, :)

  ! matrix of monthly biomass results relative
  real(RLEN), allocatable :: rel_out_monthly (:, :, :, :)
!else

!booleans for define if exists ForceNutrient and ForcePrimaryProd
  logical :: boolFN = .false.
  logical :: boolFPP = .false.

!                    *******************************
                     END MODULE DECLARATIONS_ECOPATH
!                    *******************************

