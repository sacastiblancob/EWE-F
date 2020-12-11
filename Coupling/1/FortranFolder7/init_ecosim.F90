!                    **********************
                     SUBROUTINE INIT_ECOSIM &
!                    **********************
  (TFTEL)
!***********************************************************************
! MAIN FOR CALL ECOSIM-ECOSPACE PROGRAM TROUGH JUST FEW SUBROUTINES
!***********************************************************************
!
!brief    SOLVES ECOSPACE PROGRAM
!
! Original files are part of EwE-F
! Copyright (C) 2011-2019 Middle East Technical University
! Institute of Marine Sciences (IMS-METU), Erdemli/Turkey and
! Istituto Nazionale di Oceanografia e di Geofisica Sperimentale (OGS),
! Trieste/Italy.
!
!history
! Created 31-05-2020 - Sergio Castiblanco

!#ifdef isWithBFM
!  use global_mem
!  use mem, ONLY: ppH2c,ppH1c,ppB1c,ppP1c,ppP2c,ppP3c,ppP4c, &
!       ppZ3c,ppZ4c,ppZ5c,ppZ6c,ppR6c
!  use mem, ONLY: ppH2n,ppH1n,ppB1n,ppP1n,ppP2n,ppP3n,ppP4n, &
!       ppZ3n,ppZ4n,ppZ5n,ppZ6n,ppR6n
!  use mem, ONLY: ppH2p,ppH1p,ppB1p,ppP1p,ppP2p,ppP3p,ppP4p, &
!       ppZ3p,ppZ4p,ppZ5p,ppZ6p,ppR6p
!  use mem, ONLY: ppP1l,ppP2l,ppP3l,ppP4l
!#ifdef INCLUDE_PELFE
!  use mem, ONLY: ppP1i,ppP2i,ppP3i,ppP4i
!#endif
!  use mem, ONLY: ppP1s,ppP2s,ppP3s,ppP4s
!  use mem, ONLY: iiHigherTrophicLevels
!  use statevartypesecopath, only: idxBBc, idxBBn, idxBBp, idxBBl, idxBBi, &
!       idxBBs
!  use statevartypesecosim, only: p_qpcHTL, p_qncHTL, ruHTLc, ruHTLn, ruHTLp, &
!       ruHTLl, ruHTLi, ruHTLs, tfluxc, tfluxp, tfluxn, &
!       flow2detritusR6c, flow2detritusR6p, flow2detritusR6n, &
!       multistanza_update
!#else
!  use statevartypesecopath, only: RLEN
!#endif

!  USE BIEF
!  USE DECLARATIONS_TELEMAC2D

!#ifdef _Ecospace_
  use statevartypesecospace
!#endif

  ! This is the main ecosim model. All subroutines are called from
  ! this model and state equations are calculated here.

  use statevartypesecopath, only: HDF5_fname, ep_data, ep_diet, ms_data, &
       ep_detfate, det_export_rate, flow2detritus
  use statevartypesecosim, only: GroupInfo_fname, &
       Vulnerability_fname, Forcing_fname, &
       NutrientForcing_fname, PrimaryProdForcing_fname, &
       tf, StepsPerMonth, SecondsPerMonth, NutBaseFreeProp, NutPBmax, &
       relaxeco, es_data, nvars, nstanzas, ndetritus, &
       detritus_no, BBAvg, &
       EatenOfAvg, EatenByAvg, PredAvg, LossAvg, &
       es_ms_data, imonth, FirstTime, BBeco, &
       boolFN, boolFPP, iec, j, var, &
       t0, noftsteps, step, tstep, time, integrate, &
       B, b_pred, xdot, biomeq, loss, b_out, rrate, NutTot, NutBiom, &
       NutFree, NutMin, NutFreeBase, NutBaseFreeProp, es_vul, arena, &
       frows, fcols, vrows, vcols

  use readHDF5Database

  implicit none

! INTENT VARIABLES
  REAL(RLEN),INTENT(IN)           :: TFTEL  !Here in variable is final time of Telemac

! IN  SUBROUTINE VARIABLES
  INTEGER  :: II, JJ

!  character(len = 2000):: FMT0, FMT1, FMT2
!  integer              :: vrows, vcols     ! rows & columns of vul. matrix

!#ifdef _Ecospace_
!  integer              :: lat, lon         ! coordinate variables, read from MODULE ecospace
!#endif
!  integer              :: t0               ! simulation start time, read from MODULE ecosim
!  integer              :: noftsteps, step  ! number of time steps, read from MODULE ecosim
!  real(RLEN)           :: tstep            ! time step, read from MODULE ecosim
!  real(RLEN)           :: time             ! read from MODULE ecosim

!  real(RLEN), allocatable :: B(:), b_pred(:)  ! initial biomass values       !read from MODULE ecosim
!
!  ! array of changes in biomass values calculated by derivs()
!  real(RLEN), allocatable :: xdot(:)       !read from MODULE ecosim
!!
!  ! array of changes in non-integrated biomass values calculated by derivs
!  real(RLEN), allocatable :: biomeq(:)       !read from MODULE ecosim
!
!  ! array of sinks for state variables calculated by derivs
!  real(RLEN), allocatable :: loss(:)       !read from MODULE ecosim
!
!  ! array of integrated biomass values calculated by rk4 in each time step
!  real(RLEN), allocatable :: b_out(:)       !read from MODULE ecosim
!
!  ! rate of change in state variables within the time step
!  real(RLEN), allocatable :: rrate(:)       !read from MODULE ecosim
!
!  ! (1) integrate, (0) do not integrate
!  integer, allocatable :: integrate(:)       !read from MODULE ecosim

!#ifdef _Ecospace_
!  ! matrix of biomass results absolute
!  real(RLEN), allocatable :: mat_out(:, :, :, :)      !read from MODULE ecospace
!
!  ! matrix of biomass results relative
!  real(RLEN), allocatable :: rel_out(:, :, :, :)      !read from MODULE ecospace
!
!  ! matrix of monthly biomass and catch results absolute
!  real(RLEN), allocatable :: mat_out_monthly (:, :, :, :)      !read from MODULE ecospace
!  real(RLEN), allocatable :: catch_out_monthly (:, :, :, :)      !read from MODULE ecospace
!
!  ! matrix of monthly biomass results relative
!  real(RLEN), allocatable :: rel_out_monthly (:, :, :, :)      !read from MODULE ecospace
!#else
!  ! matrix of biomass results absolute
!  real(RLEN), allocatable :: mat_out(:, :)
!
!  ! matrix of biomass results relative
!  real(RLEN), allocatable :: rel_out(:, :)
!
!  ! matrix of monthly biomass and catch results absolute
!  real(RLEN), allocatable :: mat_out_monthly (:, :)
!  real(RLEN), allocatable :: catch_out_monthly (:, :)
!
!  ! matrix of monthly biomass results relative
!  real(RLEN), allocatable :: rel_out_monthly (:, :)
!#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! get file names from namelist (nml) file
!#ifdef _Ecospace_
!  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
!       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
!       SpatialGrid_fname, SpatialDistribution_dirname, Advection_fname, &
!       ncdfout_fname, tf, SecondsPerMonth, StepsPerMonth, NutBaseFreeProp, &
!       NutPBmax, relaxeco
  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
       SpatialGrid_fname, SpatialDistribution_dirname, Advection_fname, &
       ncdfout_fname, SecondsPerMonth, StepsPerMonth, NutBaseFreeProp, &
       NutPBmax, relaxeco
!#else
!  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
!       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
!       tf, StepsPerMonth, NutBaseFreeProp, NutPBmax, relaxeco
!#endif
!  open(1010, file = "/home/aldair/Documents/EwE-F/EWE-F/Telemac-coupling-fixed/filenames.nml", status = 'OLD')
  open(1010, file = "filenames.nml", status = 'OLD')
  read(1010, nml = filenames)
  close(1010)

!!!!!!!!!! COMPUTING FINAL TIME FOR ECOSPACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tf = TFTEL/(SecondsPerMonth*12)!!

!!!!!!!!! INPUT/OUTPUT OF DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! read group info file containing Ecosim parameters
  call readEcosimScenario_io(GroupInfo_fname, nvars, nstanzas)

  ! allocate ecopath data derived type
  allocate(ep_data(nvars))

  ! allocate multistanza data derived type
  allocate(ms_data(nstanzas))

!#ifdef _Ecospace_
  ! read Spatial Grid
  call readGridFile_io(nlat, nlon, SpatialGrid_fname)

  ! read Ecospace Habitat Area Fraction (HAF) file
  call readSpatialDistribution_io(nlat, nlon, nvars, SpatialDistribution_dirname)

  ! read Advection fields
  call readAdvectionFile_io(nlat, nlon, Advection_fname)

  allocate(QperB(nlat, nlon, nvars))
  allocate(M2(nlat, nlon, nvars))
  allocate(BB_spatial(nlat, nlon, nvars))
  M2 = 0
  BB_spatial = 0
!#endif

  ! read Ecopath initial conditions file
  call readHDF5(nstanzas, HDF5_fname)

  ! read vulnerability matrix from file
  call readVulnerability_io(vrows, vcols, Vulnerability_fname)

  ! read forcing time series data from file
  call readForcingFunctions_io(Forcing_fname, nvars)

!#ifdef _ForceNutrient_
  frows = 0
  fcols = 0
  if (boolFN) then
    ! read nutrient forcing time series data from file
    call readNutrientForcingFunction_io(NutrientForcing_fname, frows)
  end if
!#endif

!#ifdef _ForcePrimaryProd_
  if (boolFPP) then
    ! read primary production forcing time series data from file
    call readPrimaryProdForcingFunction_io(PrimaryProdForcing_fname, nvars, frows, fcols)
  end if
!#endif

!#ifdef isWithBFM
!  call mapVariables (nvars)
!#endif

!!!!!!!!!! SECTION: SET SIMULATION PERIOD PARAMETERS !!!!!!!!!!!!!!!!!!!!

  t0        = 0
!#ifdef isWithBFM
!  tstep     = 1.0D0 / (12.0D0 * StepsPerMonth)
!  multistanza_update = 0.
!  imonth    = 1
!#else
  tstep     = real(1.0D0 / (12.0D0 * StepsPerMonth), RLEN)

!#endif
  WRITE(*,*) "TIIIIME", tf, tstep

  noftsteps = FLOOR(tf / tstep)
!!  noftsteps = NINT(tf / tstep)
!  IF (MOD(TFTEL/))
!  WRITE(*,*) "TSTEEEEEEP", tstep, noftsteps

ALLOCATE(ECOSUI_GRID(nlat,nlon,nvars))

!!!!!!!!!! END SECTION: SET SIMULATION PERIOD PARAMETERS !!!!!!!!!!!!!!!!

WRITE(*,*) "ECOSPACE READING DONE"

END SUBROUTINE INIT_ECOSIM
