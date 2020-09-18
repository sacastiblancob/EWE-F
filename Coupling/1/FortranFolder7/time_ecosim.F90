!                    ***************************
                     SUBROUTINE TIME_ECOSIM &
!                    ***************************
  (n,step, imonth, time, tstep, es_data, BBeco, BB_spatial, flow2detritus, &
   det_export_rate, BBAvg, LossAvg, EatenByAvg, EatenOfAvg, PredAvg, &
   es_ms_data, arena, NutFree, NutBiom, FirstTime, UpdateStanzas, &
   QperB, M2, mat_out, rel_out, COUNTER, ECOSUI, COUECO)

!
!***********************************************************************
! MAIN FOR CALL ECOSIM-ECOSPACE PROGRAM TROUGH JUST FEW SUBROUTINES
!***********************************************************************
!
!brief    COMPUTES TEMPORAL CALCS FOR ECOSPACE
!
! Original files are part of EwE-F
! Copyright (C) 2011-2019 Middle East Technical University
! Institute of Marine Sciences (IMS-METU), Erdemli/Turkey and
! Istituto Nazionale di Oceanografia e di Geofisica Sperimentale (OGS),
! Trieste/Italy.
!
!history
! Created 31-05-2020 - Sergio Castiblanco
!========================================================================

  USE BIEF
  !USE DECLARATIONS_COUPECOSPACE
  USE DECLARATIONS_TELEMAC2D, ONLY: NPOIN, ECOOUT

  use statevartypesecospace, only: nlon, nlat, grid, advection, &
       spatialhafs

!  use statevartypesecopath, only: ep_data
  use statevartypesecopath, only: RLEN, ep_data, ms_data, ep_detfate

!  use statevartypesecosim, only: BB, StepsPerMonth, UpdateStanzas, &
!       tstep, time, integrate, j, nvars, imonth, force, es_data
  use statevartypesecosim, only: ecosim_data, arena_data, &
       ecosim_multi_stanza, StepsPerMonth, &
       integrate, nvars, force, &
       nstanzas, relaxeco, noftsteps, &
       ndetritus, vrows, vcols, es_vul, &
       NutTot, NutMin, NutFreeBase, detritus_no, &
       boolFN, boolFPP

  IMPLICIT NONE

! INTENT VARIABLES
  INTEGER, INTENT(IN)   :: n, step, imonth
  REAL(RLEN), INTENT(IN) :: time, tstep
  TYPE(ecosim_data), INTENT(INOUT) :: es_data(nvars)
  REAL(RLEN), INTENT(INOUT) :: BBeco(nvars)
  REAL(RLEN), INTENT(INOUT) :: BB_spatial(nlat, nlon, nvars)
  REAL(RLEN), INTENT(INOUT), DIMENSION(ndetritus) :: flow2detritus, det_export_rate
  REAL(RLEN), INTENT(INOUT), DIMENSION(nvars) :: BBAvg, LossAvg, EatenByAvg, &
                                               EatenOfAvg, PredAvg
  TYPE(ecosim_multi_stanza), INTENT(INOUT) :: es_ms_data(nstanzas)
  TYPE(arena_data), INTENT(INOUT) :: arena
  REAL(RLEN), INTENT(INOUT)      :: NutFree, NutBiom
  LOGICAL, INTENT(INOUT)         :: FirstTime, UpdateStanzas
  REAL(RLEN), INTENT(OUT), DIMENSION(nlat, nlon, nvars) :: QperB, M2
  REAL(RLEN), INTENT(INOUT) :: mat_out(nlat, nlon, noftsteps+1, nvars)
  REAL(RLEN), INTENT(OUT) :: rel_out(nlat, nlon, noftsteps+1, nvars)

  REAL(RLEN), INTENT(IN)        :: COUNTER

  TYPE(BIEF_OBJ), INTENT(INOUT) :: ECOSUI
  TYPE(BIEF_OBJ), INTENT(IN) :: COUECO
  !TYPE(BIEF_OBJ), INTENT(OUT) :: ECOOUT
  

! IN SUBROUTINE VARIABLES
  INTEGER    :: lat, lon, latlon, j, K, II
  REAL(SELECTED_REAL_KIND(15,307)) :: CCOUN


!========================================================================
!#ifdef _Ecospace_
!!!!! run model over the specified time frame
! [...]
!      [...]
  CCOUN = 1D0/COUNTER
  !WRITE(*,*) 'CCOUN', CCOUN
! COMPUTING HABITAT SUITABILITY INDEX
  DO II = 1,nvars
    CALL OS('X=CX    ', X=ECOSUI%ADR(II)%P, C=CCOUN)
  ENDDO
  !WRITE(*,*) "ECOSUIaa", ECOSUI%ADR(1)%P%R(1)

  latlon = 1

  if (n == StepsPerMonth) then
    UpdateStanzas = .true.
  else
    UpdateStanzas = .false.
  end if

  do lon = 1, nlon
    do lat = 1, nlat
      if (grid(lat, lon) == 1) then
                          ! print instantaneous results to stdout
                          !write(*, '(A6, f9.3)' ) " Time: ", time
                          !write(*, '(A12, I4)' ) " Latitude: ", lat
                          !write(*, '(A12, I4)' ) " Longitude: ", lon

        BBeco = BB_spatial(lat, lon, :)

        call calculateFishingMortalities(nvars, imonth, BBeco, &
                                         force, ep_data, es_data)

!       do n = 1, StepsPerMonth

        ! call the Runge-Kutta 4th order numeric ode solver
        call rk4(nvars, BBeco, time, tstep, integrate, lat, lon, ep_data, &
         ms_data, ep_detfate, flow2detritus, det_export_rate, es_data, &
         nstanzas, relaxeco, StepsPerMOnth, UpdateStanzas, BBAvg, &
         LossAvg, EatenByAvg, EatenOfAvg, PredAvg, es_ms_data, ndetritus, &
         vrows, vcols, imonth, es_vul, arena, NutFree, NutTot, NutBiom, &
         NutMin, NutFreeBase, detritus_no, boolFN, boolFPP, FirstTime, &
         nlon, nlat, QperB, M2)

        ! calculate geospatial dynamics
        call ecospace(nvars, nlat, nlon, lat, lon, time, BBeco, &
                      ep_data, es_data, grid, QperB, M2, &
                      BB_spatial, advection)

        ! Update BB_spatial with new biomasses in grid
        BB_spatial(lat, lon, :) = BBeco

        mat_out(lat, lon, step + 2, :) = BBeco

        ! calculate relative change
        ! with respect to initial biomasses
        do j = 1, nvars
          rel_out(lat, lon, step + 2, j) &
                 = mat_out(lat, lon, step + 2, j) &
                 / ep_data(j)%biomass * spatialhafs(lat, lon, j)
          DO K = 1,NPOIN
            IF(INT(COUECO%R(K)) .EQ. latlon) THEN
              ECOOUT%ADR(j)%P%R(K) = mat_out(lat, lon, step + 2, j)
            ENDIF
          ENDDO
        end do

        !end do
      else
        mat_out(lat, lon, step + 2, :) = BB_spatial(lat, lon, :)
      end if
      latlon = latlon + 1
    end do
  end do
  !WRITE(*,*) 'NPOIN1', ECOOUT%ADR(7)%P%R(1)
  !WRITE(*,*) 'MATOUT', mat_out(1, 1, step + 2, 6)

! TURNING BACK TO 0 ALL HABITAT SUITABILITY VALUES
  CCOUN = 0D0;
  DO II = 1,nvars
    CALL OS('X=C     ', X=ECOSUI%ADR(II)%P, C=CCOUN)
  ENDDO
  !WRITE(*,*) "ECOSUI", ECOSUI%ADR(1)%P%R(1)

END SUBROUTINE TIME_ECOSIM
