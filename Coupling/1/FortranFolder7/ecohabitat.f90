!                    **************************
                     SUBROUTINE ECOHABITAT &
!                    **************************
   (nvars, U,V,H,ECOSUI)
!***********************************************************************
! TELEMAC2D WITH ECOSPACE
!***********************************************************************
!
! Brief. Allocates all the BIEF_OBJ for coupling with ecospace
!     + created by Sergio Castiblanco - 26/06/2020
!
!***********************************************************************

USE BIEF
USE DECLARATIONS_TELEMAC2D, ONLY: NPOIN

IMPLICIT NONE

!INTENT VARIABLES
INTEGER, INTENT(IN)           :: nvars
TYPE(BIEF_OBJ), INTENT(IN)    :: U, V, H
TYPE(BIEF_OBJ), INTENT(INOUT) :: ECOSUI

!IN SUBROUTINE VARIABLES
INTEGER                          :: II,JJ
REAL(SELECTED_REAL_KIND(15,307)) :: VEL, IND

DO II = 1,nvars
  DO JJ = 1,NPOIN
!
!   SCALAR VELOCITY    
    VEL = SQRT(U%R(JJ)**2 + V%R(JJ)**2)
!
!   SUMATORY OF INDEX
    IND = SQRT( EXP(-2*(VEL - 1.0D0)**2) * EXP(-2*(H%R(JJ) - 1.0D0)**2))
    ECOSUI%ADR(II)%P%R(JJ) = ECOSUI%ADR(II)%P%R(JJ) + IND
!
  ENDDO
ENDDO
!WRITE(*,*) "IIIIIIIND", VEL





END SUBROUTINE
