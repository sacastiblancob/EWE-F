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

program ecosim

#ifdef isWithBFM
  use global_mem
  use mem, ONLY: ppH2c,ppH1c,ppB1c,ppP1c,ppP2c,ppP3c,ppP4c, &
       ppZ3c,ppZ4c,ppZ5c,ppZ6c,ppR6c
  use mem, ONLY: ppH2n,ppH1n,ppB1n,ppP1n,ppP2n,ppP3n,ppP4n, &
       ppZ3n,ppZ4n,ppZ5n,ppZ6n,ppR6n
  use mem, ONLY: ppH2p,ppH1p,ppB1p,ppP1p,ppP2p,ppP3p,ppP4p, &
       ppZ3p,ppZ4p,ppZ5p,ppZ6p,ppR6p
  use mem, ONLY: ppP1l,ppP2l,ppP3l,ppP4l
#ifdef INCLUDE_PELFE
  use mem, ONLY: ppP1i,ppP2i,ppP3i,ppP4i
#endif
  use mem, ONLY: ppP1s,ppP2s,ppP3s,ppP4s
  use mem, ONLY: iiHigherTrophicLevels
  use statevartypesecopath, only: idxBBc, idxBBn, idxBBp, idxBBl, idxBBi, &
       idxBBs
  use statevartypesecosim, only: p_qpcHTL, p_qncHTL, ruHTLc, ruHTLn, ruHTLp, &
       ruHTLl, ruHTLi, ruHTLs, tfluxc, tfluxp, tfluxn, &
       flow2detritusR6c, flow2detritusR6p, flow2detritusR6n, &
       multistanza_update
#else
  use statevartypesecopath, only: RLEN
#endif

#ifdef _Ecospace_
  use statevartypesecospace
#endif

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

#ifdef _Ecospace_
  integer              :: lat, lon         ! coordinate variables
#endif

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

#ifdef _Ecospace_
  ! matrix of biomass results absolute
  real(RLEN), allocatable :: mat_out(:, :, :, :)

  ! matrix of biomass results relative
  real(RLEN), allocatable :: rel_out(:, :, :, :)

  ! matrix of monthly biomass and catch results absolute
  real(RLEN), allocatable :: mat_out_monthly (:, :, :, :)
  real(RLEN), allocatable :: catch_out_monthly (:, :, :, :)

  ! matrix of monthly biomass results relative
  real(RLEN), allocatable :: rel_out_monthly (:, :, :, :)
#else
  ! matrix of biomass results absolute
  real(RLEN), allocatable :: mat_out(:, :)

  ! matrix of biomass results relative
  real(RLEN), allocatable :: rel_out(:, :)

  ! matrix of monthly biomass and catch results absolute
  real(RLEN), allocatable :: mat_out_monthly (:, :)
  real(RLEN), allocatable :: catch_out_monthly (:, :)

  ! matrix of monthly biomass results relative
  real(RLEN), allocatable :: rel_out_monthly (:, :)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! get file names from namelist (nml) file
#ifdef _Ecospace_
  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
       SpatialGrid_fname, SpatialDistribution_dirname, Advection_fname, &
       tf, StepsPerMonth, NutBaseFreeProp, NutPBmax, relax
#else
  namelist /filenames/ HDF5_fname, GroupInfo_fname, Vulnerability_fname, &
       Forcing_fname, NutrientForcing_fname, PrimaryProdForcing_fname, &
       tf, StepsPerMonth, NutBaseFreeProp, NutPBmax, relax
#endif
  open(1010, file = "filenames.nml", status = 'OLD')
  read(1010, nml = filenames)
  close(1010)

!!!!!!!!!! INPUT/OUTPUT OF DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! read group info file containing Ecosim parameters
  call readEcosimScenario_io ()

  ! allocate ecopath data derived type
  allocate(ep_data(nvars))

  ! allocate multistanza data derived type
  allocate(ms_data(nstanzas))

#ifdef _Ecospace_
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
#endif

  ! read Ecopath initial conditions file
  call readHDF5 (nstanzas)

  ! read vulnerability matrix from file
  call readVulnerability_io (vrows, vcols)

  ! read forcing time series data from file
  call readForcingFunctions_io ()

#ifdef _ForceNutrient_
  ! read nutrient forcing time series data from file
  call readNutrientForcingFunction_io ()
#endif

#ifdef _ForcePrimaryProd_
  ! read primary production forcing time series data from file
  call readPrimaryProdForcingFunction_io ()
#endif

#ifdef isWithBFM
  call mapVariables (nvars)
#endif

!!!!!!!!!! SECTION: SET SIMULATION PERIOD PARAMETERS !!!!!!!!!!!!!!!!!!!!

  t0        = 0
#ifdef isWithBFM
  tstep     = 1.0D0 / (12.0D0 * StepsPerMonth)
  multistanza_update = 0.
  imonth    = 1
#else
  tstep     = real(1.0D0 / (12.0D0 * StepsPerMonth), 4)
#endif
  noftsteps = nint(tf / tstep)

!!!!!!!!!! END SECTION: SET SIMULATION PERIOD PARAMETERS !!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: MAIN CALCULATIONS BEFORE RUN MODEL !!!!!!!!!!!!!!!!!!

  ndetritus = count(ep_data(:)%org_type == 0)

  allocate(det_export_rate(ndetritus))
  allocate(detritus_no(ndetritus))
  allocate(flow2detritus(ndetritus))

  j = 0
  do i = 1, nvars
      if (ep_data(i)%org_type == 0) then
          j = j + 1
          detritus_no(j) = i
      end if
  end do

  call calculateMaximumPoBRelatedValues ()

  call calculateNutrientConcentrations ()

  call removeImportFromDiet ()

  do i = 1, nvars
      if (ep_data(i)%org_type == 2) then
          es_data(i)%Ftime = 1
          es_data(i)%hden  = es_data(i)%QB_maxoQB_0 &
               / (es_data(i)%QB_maxoQB_0 + 1)
      end if
  end do

  do i = 1, nvars
      es_data(i)%CB_base = ep_data(i)%EatenBy / ep_data(i)%biomass
      if (es_data(i)%CB_base == 0) then
          es_data(i)%CB_base = 1
          es_data(i)%Ftime_max = 1
      end if
      es_data(i)%CB_last = es_data(i)%CB_base
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
  do i = 1, nvars
      if (ep_data(i)%isstanza == 1) then
          integrate(i) = -i
      else
          integrate(i) = i
      end if
  end do

  ! setpred()
  allocate(b_pred(nvars))
  b_pred = ep_data(:)%biomass
  do i = 1, nvars
      if (ep_data(i)%biomass < 1.0e-20) then
          b_pred(i) = 1.0e-20
      end if
      if (integrate(i) >= 0) then
          es_data(i)%pred = b_pred(i)
      end if
  end do

  do i = 1, nvars
      if (es_data(i)%risk_time == -999) then
          es_data(i)%risk_time = 0
      end if
  end do

!!!!! calculate risk time for consumers
  do i = 1, nvars
      if (ep_data(i)%org_type == 2) then
          es_data(i)%CB_base = ep_data(i)%EatenBy / es_data(i)%pred
          es_data(i)%CB_last = es_data(i)%CB_base
          es_data(i)%Q_main = (1 - es_data(i)%risk_time) &
               * es_data(i)%CB_base
#ifdef isWithBFM
          es_data(i)%Q_risk = es_data(i)%risk_time * es_data(i)%CB_base &
               * (ep_data(i)%EatenOf / ep_data(i)%biomass + &
               ((1 - ep_data(i)%EE) * ep_data(i)%PoB) + 0.0000000001D0)
#else
          es_data(i)%Q_risk = real(es_data(i)%risk_time * es_data(i)%CB_base &
               * (ep_data(i)%EatenOf / ep_data(i)%biomass + &
               ((1 - ep_data(i)%EE) * ep_data(i)%PoB) + 0.0000000001D0), 4)
#endif
      end if
      WRITE(*,*) 'es_pred ', es_data(i)%pred
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

#ifndef isWithBFM


!!!!!!!!!!!!!!!!!!!!!! This is the first initalisation of Ecosim !!!!!!!!!!!!!!
  imonth = 0
  !WRITE(*,*) 'INNNNNNNN'
#ifdef _Ecospace_
  lat = 0
  lon = 0
  call derivs (0., ep_data(:)%biomass, xdot, biomeq, loss, integrate, lat, lon)
#else
  call derivs (0., ep_data(:)%biomass, xdot, biomeq, loss, integrate)
#endif
  !WRITE(*,*) 'OUT'
  FirstTime = .true.

!!!!! this is the rate of sinks in each state variable in the derivs()
  allocate (rrate(nvars))
  do i = 1, nvars
      rrate(i) = abs(loss(i)) / ep_data(i)%biomass
      !WRITE(*,*) 'VAR, rrate', i,rrate(i)
  end do

!!!!! determine whether to integrate or not each state variable
!!!!! depending on the rate of change in one time step
  do i = 1, nvars
      if (rrate(i) > 24 .and. integrate(i) == i) then
          integrate(i) = 0
          ! else if (ep_data(i)%org_type == 0) then
          !     integrate(i) = 0
      else if (ep_data(i)%isstanza == 1) then
          integrate(i) = -i
      else

      end if
      !WRITE(*,*) 'integrate', integrate(i)
  end do

!!!!!!!!!!!!!!!!!!!!!! Init is done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Now calculate final rate of change to set integration flags
  call calculateMaximumPoBRelatedValues ()
  !WRITE(*,*) 'INNNNNNNN 2'
#ifdef _Ecospace_
  lat = 0
  lon = 0
  call derivs (1., ep_data(:)%biomass, xdot, biomeq, loss, integrate, lat, lon)
#else
  call derivs (1., ep_data(:)%biomass, xdot, biomeq, loss, integrate)
#endif
  !WRITE(*,*) 'OUT 2'
  FirstTime = .true.

!!!!! this is the rate of sinks in each state variable in the derivs()
  do i = 1, nvars
      rrate(i) = abs(loss(i)) / ep_data(i)%biomass
      !WRITE(*,*) 'rrate ', rrate(i)
  end do

!!!!! determine whether to integrate or not each state variable
!!!!! depending on the rate of change in one time step
  do i = 1, nvars
      if (rrate(i) > 24 .and. integrate(i) == i) then
          integrate(i) = 0
          ! else if (ep_data(i)%org_type == 0) then
          !     integrate(i) = 0
      else if (ep_data(i)%isstanza == 1) then
          integrate(i) = -i
      else

      end if
  end do

!!!!!!!!!! SECTION: RUN MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! write the initial biomass values for the output data matrices
#ifdef _Ecospace_
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
#else
  allocate(mat_out(noftsteps + 1, nvars))
  allocate(rel_out(noftsteps + 1, nvars))
  mat_out(1, :) = ep_data(:)%biomass

  do j = 1, nvars
      rel_out(1, j) = 1.0
  end do
#endif

  ! The absolute model results will be written into the below file
  open(1111, file = "Ecosim_absResults.dat", form = 'formatted', &
       status = 'unknown')
  ! The relative model results will be written into the below file
  open(2222, file = "Ecosim_relResults.dat", form = 'formatted', &
       status = 'unknown')

  write(FMT0, '( "("I4, "(A35,"",""))" )') (nvars + 1)
  write(1111, FMT0) "Time " , groupnames
  write(2222, FMT0) "Time " , groupnames
  
  write(FMT2, '( "("I4, "(f21.7,"",""))" )') (nvars + 1)

#ifdef _Ecospace_
  write(1111, FMT2) 0.0, ep_data(:)%biomass
  write(2222, FMT2) 0.0, spread(1.0, 1, nvars)
#else
  ! Write absolute and relative model results in files
  write(1111, FMT2) 0.0, mat_out(1, :)
  write(2222, FMT2) 0.0, rel_out(1, :)
#endif

#ifdef _Ecospace_
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
#else
  allocate(mat_out_monthly((12 * tf) + 1, nvars))
  allocate(catch_out_monthly((12 * tf) + 1, nvars))
  do j = 1, nvars
      catch_out_monthly(1, j) = ep_data(j)%landings + ep_data(j)%discards
  end do

  allocate(rel_out_monthly((12 * tf) + 1, nvars))
  mat_out_monthly(1, :) = ep_data(:)%biomass

  do j = 1, nvars
      rel_out_monthly(1, j) = 1.0
  end do
#endif

  ! The absolute model results will be written into the below file
  open(3333, file = "Ecosim_absResults_monthly.dat", form = 'formatted', &
       status = 'unknown')
  ! The relative model results will be written into the below file
  open(4444, file = "Ecosim_relResults_monthly.dat", form = 'formatted', &
       status = 'unknown')
  ! The absolute monthly catches will be written into the below file
  open(5555, file = "Ecosim_absCatches_monthly.dat", form = 'formatted', &
       status = 'unknown')

  write(3333, FMT0) "Time " , groupnames
  write(4444, FMT0) "Time " , groupnames
  write(5555, FMT0) "Time " , groupnames
  
#ifdef _Ecospace_
  ! Write absolute and relative model results in files
  write(3333, FMT2) 0.0, ep_data(:)%biomass
  write(4444, FMT2) 0.0, spread(1.0, 1, nvars)
  write(5555, FMT2) 0.0, ep_data(:)%landings + ep_data(:)%discards
#else
  ! Write absolute and relative model results in files
  write(3333, FMT2) 0.0, mat_out_monthly(1, :)
  write(4444, FMT2) 0.0, rel_out_monthly(1, :)
  write(5555, FMT2) 0.0, catch_out_monthly(1, :)
#endif

  time = 0
  imonth = 0
  step = 0
#endif

  allocate(BBAvg(nvars))
  allocate(LossAvg(nvars))
  allocate(EatenByAvg(nvars))
  allocate(EatenOfAvg(nvars))
  allocate(PredAvg(nvars))

  allocate(BB(nvars))
  BB = ep_data(:)%biomass

#ifdef _Ecospace_
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
#endif
  
#ifdef isWithBFM
  allocate(ruHTLc(1,nvars,iiHigherTrophicLevels))
  allocate(ruHTLn(1,nvars,iiHigherTrophicLevels))
  allocate(ruHTLp(1,nvars,iiHigherTrophicLevels))

  allocate(ruHTLl(1,nvars,iiHigherTrophicLevels))
  allocate(ruHTLi(1,nvars,iiHigherTrophicLevels))
  allocate(ruHTLs(1,nvars,iiHigherTrophicLevels))

  allocate(tfluxc(1,iiHigherTrophicLevels))
  allocate(tfluxn(1,iiHigherTrophicLevels))
  allocate(tfluxp(1,iiHigherTrophicLevels))

  allocate(p_qncHTL(1,iiHigherTrophicLevels))
  allocate(p_qpcHTL(1,iiHigherTrophicLevels))

  p_qncHTL = (16.)/(88.5*12.)
  p_qpcHTL = (1.)/(88.5*12.)

  allocate(flow2detritusR6c(1,iiHigherTrophicLevels))
  allocate(flow2detritusR6n(1,iiHigherTrophicLevels))
  allocate(flow2detritusR6p(1,iiHigherTrophicLevels))

  print *, "Initialization ended successfully."
#endif

#ifndef isWithBFM

#ifdef _Ecospace_
!!!!! run model over the specified time frame
  do i = 0, (tf - 1)  
      do m = 1, 12        

          ! Clean monthly stanza variables
          BBAvg(:)   = 0
          LossAvg(:) = 0
          EatenByAvg(:) = 0
          EatenOfAvg(:) = 0
          PredAvg(:)    = 0

          imonth  = (i * 12) + m
          do n = 1, StepsPerMonth

              if (n == StepsPerMonth) then
                  UpdateStanzas = .true.
              else
                  UpdateStanzas = .false.
              end if
            do lon = 1, nlon
                do lat = 1, nlat
                    if (grid(lat, lon) == 1) then
                        ! this is sea
                        ! print instantaneous results to stdout
                        write(*, '(A6, f9.3)' ) " Time: ", time
                        write(*, '(A12, I4)' ) " Latitude: ", lat
                        write(*, '(A12, I4)' ) " Longitude: ", lon

                        BB = BB_spatial(lat, lon, :)

                        call calculateFishingMortalities (BB)


                        ! call the Runge-Kutta 4th order numeric ode solver
                        call rk4 (BB, time, tstep, integrate, lat, lon)

                        ! calculate geospatial dynamics
                        call ecospace (time, BB, lat, lon)

                        ! Update BB_spatial with new biomasses in grid
                        BB_spatial(lat, lon, :) = BB

                        mat_out(lat, lon, step + 2, :) = BB

                        ! calculate relative change
                        ! with respect to initial biomasses
                        do j = 1, nvars
                            rel_out(lat, lon, step + 2, j) &
                                 = mat_out(lat, lon, step + 2, j) &
                                 / ep_data(j)%biomass * spatialhafs(lat, lon, j)
                        end do

                    else
                        ! this is land
                        mat_out(lat, lon, step + 2, :) = BB_spatial(lat, lon, :)
                    end if
                
                end do
            end do
          
          step = step + 1
          time = time + tstep

          end do


      end do

  end do

  print *, "Simulation ended successfully."
  print *, "Writing netCDF file..."
  call writenetCDFfile (noftsteps, mat_out)

#else

!!!!! run model over the specified time frame
  do i = 0, (tf - 1)
  !do i = 0, 0  
      do m = 1, 12
      !do m = 1, 1        

          ! Clean monthly stanza variables
          BBAvg(:)   = 0
          LossAvg(:) = 0
          EatenByAvg(:) = 0
          EatenOfAvg(:) = 0
          PredAvg(:)    = 0

          imonth  = (i * 12) + m

          call calculateFishingMortalities (BB)
          do n = 1, StepsPerMonth
              if (n == StepsPerMonth) then
                  UpdateStanzas = .true.
              else
                  UpdateStanzas = .false.
              end if

              ! call the Runge-Kutta 4th order numeric ode solver
              call rk4 (BB, time, tstep, integrate)
              ! call EULER (BB, time, tstep, integrate)
              mat_out(step + 2, :) = BB

              ! print instantaneous results to stdout
              write(*, '(A6, f9.3)' ) " Time: ", time
              write(FMT1,  '( "(A6,", I4, "(f21.7))" )') nvars
              write(*, FMT1) "Vars: ", mat_out(step + 2, :)

              ! calculate relative change
              ! with respect to initial biomasses
              do j = 1, nvars
                  rel_out(step + 2, j) = mat_out(step + 2, j) &
                       / ep_data(j)%biomass
              end do

              ! Write absolute and relative model results in files
              write (1111, FMT2) time, mat_out(step + 2, :)
              write (2222, FMT2) time, rel_out(step + 2, :)

              step = step + 1
              time = time + tstep
          end do

          mat_out_monthly(imonth + 1, :) = BB

          do var = 1, nvars
              catch_out_monthly(imonth + 1, var) = BB(var) &
                   * es_data(var)%fishmort
          end do

          ! calculate relative change with respect to initial biomasses
          do j = 1, nvars
              rel_out_monthly(imonth + 1, j) &
                   = mat_out_monthly(imonth + 1, j) / ep_data(j)%biomass
          end do

          ! Write absolute and relative model results in files
          write(FMT2, '( "("I4, "(f27.9,"",""))" )') (nvars + 1)
          write(3333, FMT2) (time - tstep), mat_out_monthly(imonth + 1, :)
          write(4444, FMT2) (time - tstep), rel_out_monthly(imonth + 1, :)
          write(5555, FMT2) (time - tstep), catch_out_monthly(imonth + 1, :)

      end do

  end do
#endif

#endif
print *, "Simulation ended successfully."

close(1111)
close(2222)
close(3333)
close(4444)
close(5555)

deallocate(rrate)
deallocate(integrate)
deallocate(b_pred)
deallocate(xdot)
deallocate(biomeq)
deallocate(loss)
deallocate(b_out)
deallocate(mat_out)
deallocate(rel_out)
deallocate(mat_out_monthly)
deallocate(catch_out_monthly)
deallocate(rel_out_monthly)

call freeMemory ()

end program
