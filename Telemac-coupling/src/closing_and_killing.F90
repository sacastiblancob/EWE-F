!                    ******************************
                     SUBROUTINE CLOSING_AND_KILLING
!                    ******************************
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

  use statevartypesecopath

  use statevartypesecosim

  IMPLICIT NONE

!  close(1111)
!  close(2222)
!  close(3333)
!  close(4444)
!  close(5555)

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



END SUBROUTINE CLOSING_AND_KILLING
