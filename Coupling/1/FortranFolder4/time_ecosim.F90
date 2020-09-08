!                    ***************************
                     SUBROUTINE TIME_ECOSIM(n, step)
!                    ***************************
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

  use statevartypesecospace

  use statevartypesecopath, only: ep_data

  use statevartypesecosim, only: BB, StepsPerMonth, UpdateStanzas, &
       tstep, time, integrate, j, nvars


  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: n, step

!========================================================================
!#ifdef _Ecospace_
!!!!! run model over the specified time frame
! [...]
!      [...]
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

                      BB = BB_spatial(lat, lon, :)

                      call calculateFishingMortalities (BB)

!                     do n = 1, StepsPerMonth

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

                      !end do
                  else
                      mat_out(lat, lon, step + 2, :) = BB_spatial(lat, lon, :)
                  end if

              end do
          end do
          WRITE(*,*) "MAT_OUT", mat_out(2, 2, step + 2, :)

!      [...]
! [...]
!#else
!
!!!!!! run model over the specified time frame
!  [...]
!      [...]
!
!          ! Clean monthly stanza variables
!          BBAvg(:)   = 0
!          LossAvg(:) = 0
!          EatenByAvg(:) = 0
!          EatenOfAvg(:) = 0
!          PredAvg(:)    = 0
!
!          imonth  = (i * 12) + m
!
!          call calculateFishingMortalities (BB)
!          do n = 1, StepsPerMonth
!              if (n == StepsPerMonth) then
!                  UpdateStanzas = .true.
!              else
!                  UpdateStanzas = .false.
!              end if
!
!              ! call the Runge-Kutta 4th order numeric ode solver
!              call rk4 (BB, time, tstep, integrate)
!              mat_out(step + 2, :) = BB
!
!              ! print instantaneous results to stdout
!              write(*, '(A6, f9.3)' ) " Time: ", time
!              write(FMT1,  '( "(A6,", I4, "(f21.7))" )') nvars
!              write(*, FMT1) "Vars: ", mat_out(step + 2, :)
!
!              ! calculate relative change
!              ! with respect to initial biomasses
!              do j = 1, nvars
!                  rel_out(step + 2, j) = mat_out(step + 2, j) &
!                       / ep_data(j)%biomass
!              end do
!
!              ! Write absolute and relative model results in files
!              write (1111, FMT2) time, mat_out(step + 2, :)
!              write (2222, FMT2) time, rel_out(step + 2, :)
!
!              step = step + 1
!              time = time + tstep
!          end do
!
!          mat_out_monthly(imonth + 1, :) = BB
!
!          do var = 1, nvars
!              catch_out_monthly(imonth + 1, var) = BB(var) &
!                   * es_data(var)%fishmort
!          end do
!
!          ! calculate relative change with respect to initial biomasses
!          do j = 1, nvars
!              rel_out_monthly(imonth + 1, j) &
!                   = mat_out_monthly(imonth + 1, j) / ep_data(j)%biomass
!          end do
!
!          ! Write absolute and relative model results in files
!          write(FMT2, '( "("I4, "(f27.9,"",""))" )') (nvars + 1)
!          write(3333, FMT2) (time - tstep), mat_out_monthly(imonth + 1, :)
!          write(4444, FMT2) (time - tstep), rel_out_monthly(imonth + 1, :)
!          write(5555, FMT2) (time - tstep), catch_out_monthly(imonth + 1, :)
!
!      [...]
!
!  [...]
!#endif




END SUBROUTINE TIME_ECOSIM
