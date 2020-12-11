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
REAL(SELECTED_REAL_KIND(15,307)) :: VEL, INDH, INDV

DO II = 1,nvars
  DO JJ = 1,NPOIN
		IF ((H%R(JJ).LT.0.7D0) .AND. (H%R(JJ).GT.0.0D0)) THEN
		  INDH=0.00000236D0*((100*H%R(JJ))**3)-0.000578D0 &
          *((100*H%R(JJ))**2)+0.024D0*(100*H%R(JJ))+0.442D0
		ELSEIF (H%R(JJ).LE.0.0D0) THEN
		  INDH=0.0D0
		ELSE
		  INDH=0.0993D0
		ENDIF
!
!   SCALAR VELOCITY    
    VEL = SQRT(U%R(JJ)**2 + V%R(JJ)**2)
!
		IF ((VEL.LT.1.0D0) .AND. (VEL.GT.0.2D0)) THEN
      INDV=-0.023D0*(VEL**3)+0.508D0*(VEL**2) &
         +0.3755D0*(VEL)+0.1213D0
    ELSEIF (VEL.LE.0.2D0) THEN
      INDV=0.2165D0
    ELSE
      INDV=0.9818D0
    ENDIF
!   SUMATORY OF INDEX
!    IND = SQRT( EXP(-2*(VEL - 1.0D0)**2) * EXP(-2*(H%R(JJ) - 1.0D0)**2))
    ECOSUI%ADR(II)%P%R(JJ) = ECOSUI%ADR(II)%P%R(JJ) + SQRT(INDH*INDV)
    INDH = 0D0
    INDV = 0D0
!
  ENDDO
ENDDO
!WRITE(*,*) "IIIIIIIND", VEL





END SUBROUTINE
