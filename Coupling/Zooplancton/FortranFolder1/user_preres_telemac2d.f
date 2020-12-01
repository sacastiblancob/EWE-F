!                    ********************************
                     SUBROUTINE USER_PRERES_TELEMAC2D
!                    ********************************
!
!***********************************************************************
! TELEMAC2D
!***********************************************************************
!
!brief    PREPARES THE USER VARIABLES WHICH WILL BE WRITTEN TO
!+                THE RESULTS FILE OR TO THE LISTING.
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D, EX_USER_PRERES_TELEMAC2D
     &                         => USER_PRERES_TELEMAC2D

!
      USE DECLARATIONS_SPECIAL
      USE INTERFACE_PARALLEL, ONLY : P_DMAX,P_DMIN
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER N
!	  DOUBLE PRECISION TMP1, TMP2
!
!=======================================================================
! COMPUTES SOME INDEX 0
!=======================================================================
!
      IF((LEO.AND.SORLEO(23)).OR.(IMP.AND.SORIMP(23))) THEN
	    DO N=1,NPOIN
!		  IF ((H%R(N).LT.0.7D0) .AND. (H%R(N).GT.0.0D0)) THEN
!!		    TMP1=-0.000578D0*((100*H%R(N))**2)
!!			TMP2=0.024D0*(100*H%R(N))+0.442D0+100.0D0
!		    PRIVE%ADR(1)%P%R(N)=0.00000236D0*((100*H%R(N))**3)-0.000578D0
!     &      *((100*H%R(N))**2)+0.024D0*(100*H%R(N))+0.442D0
!		  ELSEIF (H%R(N).LE.0.0D0) THEN
!		    PRIVE%ADR(1)%P%R(N)=0.0D0
!		  ELSE
!		    PRIVE%ADR(1)%P%R(N)=0.0993D0
!		  ENDIF
            PRIVE%ADR(1)%P%R(N)=3.0D0
        ENDDO
      ENDIF
!
!=======================================================================
! COMPUTES SOME INDEX 1
!=======================================================================
!
      IF((LEO.AND.SORLEO(24)).OR.(IMP.AND.SORIMP(24))) THEN
		DO N=1,NPOIN
              !PRIVE%ADR(2)%P%R(N)=T6%R(N)+2.0D0
              PRIVE%ADR(2)%P%R(N)=2.0D0
        ENDDO
      ENDIF
!
!=======================================================================
! COMPUTES SOME INDEX 2
!=======================================================================
!
      IF((LEO.AND.SORLEO(24)).OR.(IMP.AND.SORIMP(24))) THEN
		DO N=1,NPOIN
              !PRIVE%ADR(3)%P%R(N)=SQRT(U%R(N)**2+V%R(N)**2)+3.0D0
              PRIVE%ADR(3)%P%R(N)=1.0D0
        ENDDO
      ENDIF
!
!=======================================================================
! COMPUTES SOME INDEX 3
!=======================================================================
!
      IF((LEO.AND.SORLEO(24)).OR.(IMP.AND.SORIMP(24))) THEN
        DO N=1,NPOIN
!          PRIVE%ADR(4)%P%R(N)=SQRT(U%R(N)**2+V%R(N)**2)+4.0D0
          PRIVE%ADR(4)%P%R(N)=0.5D0
      ENDDO
      ENDIF
!
!=======================================================================
! COMPUTES SOME INDEX 4 (test!!!!!!!)
!=======================================================================
!
!     PRINT *, MAXVAL(IDX%R)
!      GRUF1 = ECOUT%ADR(6)
      DO N=1,NPOIN
!	    GRUF1%R(N)=ECOOUT%ADR(6)%P%R(N)
          GRUF1%R(N) = 1.0D0
      ENDDO
!
!=======================================================================
!
      RETURN
      END
