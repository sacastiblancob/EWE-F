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

  use statevartypesecosim, only: tf, iec, m, noftsteps

  use statevartypesecospace, only: mat_out, ncdfout_fname

  IMPLICIT NONE

!================================================================================

  CALL INIT_ECOSIM()

  WRITE(*,*) "FINAL TIME", tf

  do iec = 0, (tf - 1)
    do m = 1, 12
      CALL TIME_ECOSIM(iec,m)
    end do
  end do

  print *, "Writing netCDF file in ", ncdfout_fname
  call writenetCDFfile (noftsteps, mat_out, ncdfout_fname)

  print *, "Simulation ended successfully."

  CALL CLOSING_AND_KILLING()


END PROGRAM MAIN
