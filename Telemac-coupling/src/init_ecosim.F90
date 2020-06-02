!                    **********************
                     SUBROUTINE INIT_ECOSIM
!                    **********************
!
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

!#ifdef _Ecospace_
  use statevartypesecospace
!#endif

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
       es_ms_data, imonth, FirstTime, UpdateStanzas, BB, groupnames, &
       boolFN, boolFPP, iec, j, var, &
       t0, noftsteps, step, tstep, time, integrate, &
       B, b_pred, xdot, biomeq, loss, b_out, rrate
  use readHDF5Database

  implicit none

  character(len = 2000):: FMT0, FMT1, FMT2
  integer              :: vrows, vcols     ! rows & columns of vul. matrix

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
  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
       SpatialGrid_fname, SpatialDistribution_dirname, Advection_fname, &
       ncdfout_fname, tf, StepsPerMonth, NutBaseFreeProp, NutPBmax, relax
!#else
!  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
!       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
!       tf, StepsPerMonth, NutBaseFreeProp, NutPBmax, relax
!#endif
  open(1010, file = "/home/aldair/Documents/EwE-F/git/Telemac-coupling/filenames.nml", status = 'OLD')
  read(1010, nml = filenames)
  close(1010)

!!!!!!!!!! INPUT/OUTPUT OF DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! read group info file containing Ecosim parameters
  call readEcosimScenario_io()

  ! allocate ecopath data derived type
  allocate(ep_data(nvars))

  ! allocate multistanza data derived type
  allocate(ms_data(nstanzas))

!#ifdef _Ecospace_
  ! read Spatial Grid
  call readGridFile_io ()

  ! read Ecospace Habitat Area Fraction (HAF) file
  call readSpatialDistribution_io ()

  ! read Advection fields
  call readAdvectionFile_io()

  allocate(QperB(nlat, nlon, nvars))
  allocate(M2(nlat, nlon, nvars))
  allocate(BB_spatial(nlat, nlon, nvars))
  M2 = 0
  BB_spatial = 0
!#endif

  ! read Ecopath initial conditions file
  call readHDF5 (nstanzas)

  ! read vulnerability matrix from file
  call readVulnerability_io (vrows, vcols)

  ! read forcing time series data from file
  call readForcingFunctions_io ()

!#ifdef _ForceNutrient_
  if (boolFN) then
    ! read nutrient forcing time series data from file
    call readNutrientForcingFunction_io ()
  end if
!#endif

!#ifdef _ForcePrimaryProd_
  if (boolFPP) then
    ! read primary production forcing time series data from file
    call readPrimaryProdForcingFunction_io ()
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
  tstep     = real(1.0D0 / (12.0D0 * StepsPerMonth), 4)
!#endif
  noftsteps = nint(tf / tstep)

!!!!!!!!!! END SECTION: SET SIMULATION PERIOD PARAMETERS !!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: MAIN CALCULATIONS BEFORE RUN MODEL !!!!!!!!!!!!!!!!!!

  ndetritus = count(ep_data(:)%org_type == 0)

  allocate(det_export_rate(ndetritus))
  allocate(detritus_no(ndetritus))
  allocate(flow2detritus(ndetritus))

  j = 0
  do iec = 1, nvars
      if (ep_data(iec)%org_type == 0) then
          j = j + 1
          detritus_no(j) = iec
      end if
  end do

  call calculateMaximumPoBRelatedValues ()

  call calculateNutrientConcentrations ()

  call removeImportFromDiet ()

  do iec = 1, nvars
      if (ep_data(iec)%org_type == 2) then
          es_data(iec)%Ftime = 1
          es_data(iec)%hden  = es_data(iec)%QB_maxoQB_0 &
               / (es_data(iec)%QB_maxoQB_0 + 1)
      end if
!      write(*,*) es_data(iec)%QB_maxoQB_0
!      write(*,*) es_data(iec)%hden
  end do

  do iec = 1, nvars
      es_data(iec)%CB_base = ep_data(iec)%EatenBy / ep_data(iec)%biomass
      if (es_data(iec)%CB_base == 0) then
          es_data(iec)%CB_base = 1
          es_data(iec)%Ftime_max = 1
      end if
      es_data(iec)%CB_last = es_data(iec)%CB_base
!      write(*,*) i, es_data(iec)%CB_base
!      write(*,*) es_data(iec)%Ftime_max
!      write(*,*) es_data(iec)%CB_last
!      write(*,*) ep_data(iec)%EatenBy
!      write(*,*) ep_data(iec)%biomass
!      write(*,*)
  end do

!!!!! initialize stanza parameters if any
  allocate(B(nvars))
  B = 0
  if (nstanzas > 0) then
      allocate (es_ms_data(nstanzas))
      call initialiseSplitGroups (B)
  end if
  deallocate(B)

  allocate(integrate(nvars))
  do iec = 1, nvars
      if (ep_data(iec)%isstanza == 1) then
          integrate(iec) = -iec
      else
          integrate(iec) = iec
      end if
!      write(*,*) i,ep_data(iec)%isstanza
  end do

  ! setpred()
  allocate(b_pred(nvars))
  b_pred = ep_data(:)%biomass
  do iec = 1, nvars
!      write(*,*) i,es_data(iec)%pred
      if (ep_data(iec)%biomass < 1.0e-20) then
          b_pred(iec) = 1.0e-20
      end if
      if (integrate(iec) >= 0) then
          es_data(iec)%pred = b_pred(iec)
      end if
!      write(*,*) i,es_data(iec)%pred
  end do

  do iec = 1, nvars
!      write(*,*) es_data(iec)%risk_time
      if (es_data(iec)%risk_time == -999) then
          es_data(iec)%risk_time = 0
      end if
  end do

!!!!! calculate risk time for consumers
  do iec = 1, nvars
      if (ep_data(iec)%org_type == 2) then
          es_data(iec)%CB_base = ep_data(iec)%EatenBy / es_data(iec)%pred
          es_data(iec)%CB_last = es_data(iec)%CB_base
          es_data(iec)%Q_main = (1 - es_data(iec)%risk_time) &
               * es_data(iec)%CB_base
!#ifdef isWithBFM
!          es_data(iec)%Q_risk = es_data(iec)%risk_time * es_data(iec)%CB_base &
!               * (ep_data(iec)%EatenOf / ep_data(iec)%biomass + &
!               ((1 - ep_data(iec)%EE) * ep_data(iec)%PoB) + 0.0000000001D0)
!#else
          es_data(iec)%Q_risk = real(es_data(iec)%risk_time * es_data(iec)%CB_base &
               * (ep_data(iec)%EatenOf / ep_data(iec)%biomass + &
               ((1 - ep_data(iec)%EE) * ep_data(iec)%PoB) + 0.0000000001D0), 4)
!#endif
      end if
  end do

  call calculateLotkaVolterraEffectiveSearchRates (vrows, vcols)

!!!!! initial switching parameters (InitRelaSwitch)
  call initialiseRelativeSwitchingParameters (vrows, vcols)

!!!!! set arena vulnerability and search rates
  call setArenaVulnerabilityandSearchRates (vrows, vcols)

!!!!! also allocate output variables from derivs()
  allocate(xdot(nvars))
  allocate(biomeq(nvars))
  allocate(loss(nvars))
  allocate(b_out(nvars))

!!!!!! prepare for run model by calculating initial values
  ! for preparing initial values, year is set to null
  ! to disregard forcing data

  do j = 1, nvars
      es_data(j)%fishmort = ((ep_data(j)%landings + ep_data(j)%discards) &
           / ep_data(j)%biomass)
  end do

  FirstTime = .true.

!#ifndef isWithBFM


!!!!!!!!!!!!!!!!!!!!!! This is the first initalisation of Ecosim !!!!!!!!!!!!!!
  imonth = 0
!#ifdef _Ecospace_
  lat = 0
  lon = 0
!  call derivs (0., ep_data(:)%biomass, xdot, biomeq, loss, integrate, lat, lon) !Original, 0. was taken weird
  call derivs (0D0, ep_data(:)%biomass, xdot, biomeq, loss, integrate, lat, lon)
!#else
!!  call derivs (1., ep_data(:)%biomass, xdot, biomeq, loss, integrate) !Original, 0. was taken weird
!  call derivs (0D0, ep_data(:)%biomass, xdot, biomeq, loss, integrate)
!#endif
  FirstTime = .true.

!!!!! this is the rate of sinks in each state variable in the derivs()
  allocate (rrate(nvars))
  do iec = 1, nvars
      rrate(iec) = abs(loss(iec)) / ep_data(iec)%biomass
  end do

!!!!! determine whether to integrate or not each state variable
!!!!! depending on the rate of change in one time step
  do iec = 1, nvars
      if (rrate(iec) > 24 .and. integrate(iec) == iec) then
          integrate(iec) = 0
          ! else if (ep_data(iec)%org_type == 0) then
          !     integrate(iec) = 0
      else if (ep_data(iec)%isstanza == 1) then
          integrate(iec) = -iec
      else

      end if
  end do

!!!!!!!!!!!!!!!!!!!!!! Init is done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Now calculate final rate of change to set integration flags
  call calculateMaximumPoBRelatedValues ()
!#ifdef _Ecospace_
  lat = 0
  lon = 0
!  call derivs (1., ep_data(:)%biomass, xdot, biomeq, loss, integrate, lat, lon) !Original, 1. was taken weird
  call derivs (1D0, ep_data(:)%biomass, xdot, biomeq, loss, integrate, lat, lon)
!#else
!!  call derivs (1., ep_data(:)%biomass, xdot, biomeq, loss, integrate) !Original, 1. was taken weird
!  call derivs (1D0, ep_data(:)%biomass, xdot, biomeq, loss, integrate)
!#endif

  FirstTime = .true.

!!!!! this is the rate of sinks in each state variable in the derivs()
  do iec = 1, nvars
      rrate(iec) = abs(loss(iec)) / ep_data(iec)%biomass
  end do

!!!!! determine whether to integrate or not each state variable
!!!!! depending on the rate of change in one time step
  do iec = 1, nvars
      if (rrate(iec) > 24 .and. integrate(iec) == iec) then
          integrate(iec) = 0
          ! else if (ep_data(iec)%org_type == 0) then
          !     integrate(iec) = 0
      else if (ep_data(iec)%isstanza == 1) then
          integrate(iec) = -iec
      else

      end if
  end do


!!!!!!!!!! SECTION: RUN MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! write the initial biomass values for the output data matrices
!#ifdef _Ecospace_
  allocate(mat_out(nlat, nlon, noftsteps + 1, nvars))
  allocate(rel_out(nlat, nlon, noftsteps + 1, nvars))
  do lat = 1, nlat
      do lon = 1, nlon
          do var = 1, nvars
              mat_out(lat, lon, 1, var) = ep_data(var)%biomass &
                   * spatialhafs(lat, lon, var)
          end do
      end do
  end do

  do j = 1, nvars
      rel_out(:, :, 1, j) = 1.0
  end do
!#else
!  allocate(mat_out(noftsteps + 1, nvars))
!  allocate(rel_out(noftsteps + 1, nvars))
!  mat_out(1, :) = ep_data(:)%biomass

!  do j = 1, nvars
!      rel_out(1, j) = 1.0
!  end do
!#endif

  ! The absolute model results will be written into the below file
!  open(1111, file = "Ecosim_absResults.dat", form = 'formatted', &
!       status = 'unknown')
  ! The relative model results will be written into the below file
!  open(2222, file = "Ecosim_relResults.dat", form = 'formatted', &
!       status = 'unknown')

!  write(FMT0, '( "("I4, "(A35,"",""))" )') (nvars + 1)
!  write(1111, FMT0) "Time " , groupnames
!  write(2222, FMT0) "Time " , groupnames

!  write(FMT2, '( "("I4, "(f21.7,"",""))" )') (nvars + 1)

!#ifdef _Ecospace_
!  write(1111, FMT2) 0.0, ep_data(:)%biomass
!  write(2222, FMT2) 0.0, spread(1.0, 1, nvars)
!#else
!  ! Write absolute and relative model results in files
!  write(1111, FMT2) 0.0, mat_out(1, :)
!  write(2222, FMT2) 0.0, rel_out(1, :)
!#endif

!#ifdef _Ecospace_
  allocate(mat_out_monthly(nlat, nlon, (12 * tf) + 1, nvars))
  do lat = 1, nlat
      do lon = 1, nlon
          do var = 1, nvars
              mat_out_monthly(lat, lon, 1, var) = ep_data(var)%biomass &
                   * spatialhafs(lat, lon, var)
          end do
      end do
  end do

  allocate(catch_out_monthly(nlat, nlon, (12 * tf) + 1, nvars))
  do lat = 1, nlat
      do lon = 1, nlon
          do j = 1, nvars
              catch_out_monthly(lat, lon, 1, j) = ep_data(j)%landings &
                   * spatialhafs(lat, lon, j) &
                   + ep_data(j)%discards * spatialhafs(lat, lon, j)
          end do
      end do
  end do

  allocate(rel_out_monthly(nlat, nlon, (12 * tf) + 1, nvars))
  do j = 1, nvars
      rel_out_monthly(:, :, 1, j) = 1.0
  end do
!#else
!  allocate(mat_out_monthly((12 * tf) + 1, nvars))
!  allocate(catch_out_monthly((12 * tf) + 1, nvars))
!  do j = 1, nvars
!      catch_out_monthly(1, j) = ep_data(j)%landings + ep_data(j)%discards
!  end do
!
!  allocate(rel_out_monthly((12 * tf) + 1, nvars))
!  mat_out_monthly(1, :) = ep_data(:)%biomass
!
!  do j = 1, nvars
!      rel_out_monthly(1, j) = 1.0
!  end do
!#endif

  ! The absolute model results will be written into the below file
!  open(3333, file = "Ecosim_absResults_monthly.dat", form = 'formatted', &
!       status = 'unknown')
  ! The relative model results will be written into the below file
!  open(4444, file = "Ecosim_relResults_monthly.dat", form = 'formatted', &
!       status = 'unknown')
  ! The absolute monthly catches will be written into the below file
!  open(5555, file = "Ecosim_absCatches_monthly.dat", form = 'formatted', &
!       status = 'unknown')

!  write(3333, FMT0) "Time " , groupnames
!  write(4444, FMT0) "Time " , groupnames
!  write(5555, FMT0) "Time " , groupnames

!#ifdef _Ecospace_
  ! Write absolute and relative model results in files
!  write(3333, FMT2) 0.0, ep_data(:)%biomass
!  write(4444, FMT2) 0.0, spread(1.0, 1, nvars)
!  write(5555, FMT2) 0.0, ep_data(:)%landings + ep_data(:)%discards
!#else
!  ! Write absolute and relative model results in files
!  write(3333, FMT2) 0.0, mat_out_monthly(1, :)
!  write(4444, FMT2) 0.0, rel_out_monthly(1, :)
!  write(5555, FMT2) 0.0, catch_out_monthly(1, :)
!#endif

  time = 0
  imonth = 0
  step = 0
!  WRITE(*,*) 'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
!#endif

!  write(*,*) 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'


  allocate(BBAvg(nvars))
  allocate(LossAvg(nvars))
  allocate(EatenByAvg(nvars))
  allocate(EatenOfAvg(nvars))
  allocate(PredAvg(nvars))

  allocate(BB(nvars))
  BB = ep_data(:)%biomass

!#ifdef _Ecospace_
  do lon = 1, nlon
      do lat = 1, nlat
          do var = 1, nvars
              if (grid(lat, lon) == 1) then
                  BB_spatial(lat, lon, var) = BB(var) * spatialhafs(lat, lon, var)
              else
                  BB_spatial(lat, lon, var) = -999
              end if
          end do
      end do
  end do
!#endif

!#ifdef isWithBFM
!  allocate(ruHTLc(1,nvars,iiHigherTrophicLevels))
!  allocate(ruHTLn(1,nvars,iiHigherTrophicLevels))
!  allocate(ruHTLp(1,nvars,iiHigherTrophicLevels))
!
!  allocate(ruHTLl(1,nvars,iiHigherTrophicLevels))
!  allocate(ruHTLi(1,nvars,iiHigherTrophicL!evels))
!  allocate(ruHTLs(1,nvars,iiHigherTrophicLevels))
!
!  allocate(tfluxc(1,iiHigherTrophicLevels))
!  allocate(tfluxn(1,iiHigherTrophicLevels))
!  allocate(tfluxp(1,iiHigherTrophicLevels))
!
!  allocate(p_qncHTL(1,iiHigherTrophicLevels))
!  allocate(p_qpcHTL(1,iiHigherTrophicLevels))
!
!  p_qncHTL = (16.)/(88.5*12.)
!  p_qpcHTL = (1.)/(88.5*12.)
!
!  allocate(flow2detritusR6c(1,iiHigherTrophicLevels))
!  allocate(flow2detritusR6n(1,iiHigherTrophicLevels))
!  allocate(flow2detritusR6p(1,iiHigherTrophicLevels))
!
!  print *, "Initialization ended successfully."
!#endif

!  WRITE(*,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'

!#ifndef isWithBFM

!  WRITE(*,*) 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD'

WRITE(*,*) "FINAL TIME", tf
WRITE(*,*) "ECOSPACE INITIALISATION DONE"



END SUBROUTINE INIT_ECOSIM
