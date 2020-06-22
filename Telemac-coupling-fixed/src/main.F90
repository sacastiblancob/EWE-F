!                    ************
                     PROGRAM MAIN
!                    ************
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

  use statevartypesecopath, only: RLEN, flow2detritus, det_export_rate

  use statevartypesecosim, only: BBAvg, LossAvg, EatenByAvg, EatenOfAvg, &
       PredAvg, imonth, tstep, step, time, tf, StepsPerMonth, SecondsPerMonth, &
       iec, n, noftsteps, es_data, BB, BBAvg, LossAvg, EatenByAvg, EatenOfAvg, &
       PredAvg, es_ms_data, arena, NutFree, NutBiom, FirstTime, UpdateStanzas, &
       noftsteps

  use statevartypesecospace

  IMPLICIT NONE

  !Variables that emulate variables from telemac
  REAL(RLEN)      :: TFTEL   !FINAL TIME IN TELEMAC
  REAL(RLEN)      :: TSTEL    !TIME STEP IN TELEMAC
  REAL(RLEN)      :: TIMETEL  !TIME EXECUTION IN TELEMAC
  INTEGER         :: NTSTEL   !NUMBER OF TIME STEPS IN TELEMAC

!================================================================================

  TFTEL = 2401.0D0
  TSTEL = 1.0
  TIMETEL = 0
  NTSTEL = INT(TFTEL/TSTEL)

  CALL INIT_ECOSIM(TFTEL)

  WRITE(*,*) "FINAL TIME", tf
  imonth = 1
!  step = 0
  n = 1

  WRITE(*,*) "NNNNNOFTSTEPS", noftsteps

  do iec = 0,NTSTEL
!  do iec = 0, (tf - 1)
!    do m = 1, 12
   IF(TIMETEL.GT.(SecondsPerMonth*imonth)) THEN

       WRITE(*,*) step, imonth, n

!       CALL TIME_ECOSIM(n,step)
!        CALL TIME_ECOSIM(n,step, imonth, time, tstep, es_data, BB, BB_spatial,
!     &                   flow2detritus, det_export_rate, BBAvg, LossAvg,
!     &                   EatenByAvg, EatenOfAvg, PredAvg, es_ms_data, arena,
!     &                   NutFree, NutBiom, FirstTime, UpdateStanzas, QperB,
!     &                   M2, mat_out, rel_out)
        CALL TIME_ECOSIM(n,step, imonth, time, tstep, es_data, BB, BB_spatial, &
                        flow2detritus, det_export_rate, BBAvg, LossAvg, &
                        EatenByAvg, EatenOfAvg, PredAvg, es_ms_data, arena, &
                        NutFree, NutBiom, FirstTime, UpdateStanzas, QperB, &
                        M2, mat_out, rel_out)
       step = step + 1
       time = time + tstep
       n = n + 1

       ! Clean monthly stanza variables
       BBAvg(:)   = 0
       LossAvg(:) = 0
       EatenByAvg(:) = 0
       EatenOfAvg(:) = 0
       PredAvg(:)    = 0

       imonth = imonth + 1
       n = 1
       write(*,*) tf, imonth, time, step

   ENDIF

   IF(TIMETEL.GT.(SecondsPerMonth*(imonth-1) &
        + n*(SecondsPerMonth/StepsPerMonth))) THEN
        WRITE(*,*) step, imonth, n
!       imonth  = (iec * 12) + m
!      do n = 1, StepsPerMonth
!         write(*,*) iec, m, n, imonth, time

!       CALL TIME_ECOSIM(n,step)
!        CALL TIME_ECOSIM(n,step, imonth, time, tstep, es_data, BB, BB_spatial,
!     &                   flow2detritus, det_export_rate, BBAvg, LossAvg,
!     &                   EatenByAvg, EatenOfAvg, PredAvg, es_ms_data, arena,
!     &                   NutFree, NutBiom, FirstTime, UpdateStanzas, QperB,
!     &                   M2, mat_out, rel_out)
        CALL TIME_ECOSIM(n,step, imonth, time, tstep, es_data, BB, BB_spatial, &
                        flow2detritus, det_export_rate, BBAvg, LossAvg, &
                        EatenByAvg, EatenOfAvg, PredAvg, es_ms_data, arena, &
                        NutFree, NutBiom, FirstTime, UpdateStanzas, QperB, &
                        M2, mat_out, rel_out)

         step = step + 1
         time = time + tstep
         n = n + 1

!       end do
   ENDIF
!    end do
    TIMETEL = TIMETEL+TSTEL
  end do

  print *, "Writing netCDF file in ", ncdfout_fname
  call writenetCDFfile (noftsteps, mat_out, ncdfout_fname)

  print *, "Simulation ended successfully."

!  CALL CLOSING_AND_KILLING()


END PROGRAM MAIN
