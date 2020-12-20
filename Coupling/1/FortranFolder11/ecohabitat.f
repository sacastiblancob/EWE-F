!                    **************************
                     SUBROUTINE ECOHABITAT
!                    **************************
     &(nvars, U,V,H,ECOSUI,ECGRXN,ECGRYN,UNSV2D,MESH,MSK,MASKEL,S)
!***********************************************************************
! TELEMAC2D WITH ECOSPACE
!***********************************************************************
!
! Brief. Allocates all the BIEF_OBJ for coupling with ecospace
!     + created by Sergio Castiblanco - 26/06/2020
!
!***********************************************************************
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| nvars          |-->| NUMBER OF TROPHIC GROUPS
!| U              |-->| VELOCITY X FIELD
!| V              |-->| VELOCITY Y FIELD
!| H              |-->| WATER DEPTH
!| ECOSUI         |<->| HABITAT SUITABILITY FIELD
!| ECGRAX         |<->| GRADIENT OF HABITAT SUITABILITY FIELD IN X
!| ECGRAY         |<->| GRADIENT OF HABITAT SUITABILITY FIELD IN Y
!| UNSV2D         |-->| INVERSE OF INTEGRALS OF FEM TEST FUNCTIONS
!| MESH           |-->| MESH STRUCTURE
!| MSK            |-->| IF YES, THERE IS MASKED ELEMENTS.
!| MASKEL         |-->| MASKING OF ELEMENTS
!| S              |-->| VOID STRUCTURE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      USE BIEF
      USE DECLARATIONS_TELEMAC2D, ONLY: NPOIN

      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!INTENT VARIABLES
      INTEGER, INTENT(IN)           :: nvars
      LOGICAL, INTENT(IN)           :: MSK
!
!  STRUCTURE OF VECTORS
!
      TYPE(BIEF_OBJ), INTENT(IN)    :: U, V, H, UNSV2D, MASKEL
      TYPE(BIEF_OBJ), INTENT(INOUT) :: ECGRXN, ECGRYN
!     DUMMY STRUCTURE
      TYPE(BIEF_OBJ), INTENT(IN)    :: S
!
!  STRUCTURE OF BLOCKS
!
      TYPE(BIEF_OBJ), INTENT(INOUT) :: ECOSUI
!
!  STRUCTURE OF MESH
!
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!IN SUBROUTINE VARIABLES
      INTEGER                          :: II,JJ,IELM
      REAL(SELECTED_REAL_KIND(15,307)) :: VEL, INDH, INDV
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
! COMPUTING HABITAT SUITABILITY, AND STORED IN ECOSUI BLOCK
!
      DO II = 1,nvars
       DO JJ = 1,NPOIN
		IF ((H%R(JJ).LT.0.7D0) .AND. (H%R(JJ).GT.0.0D0)) THEN
		 INDH=0.00000236D0*((100*H%R(JJ))**3)-0.000578D0
     &    *((100*H%R(JJ))**2)+0.024D0*(100*H%R(JJ))+0.442D0
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
         INDV=-0.023D0*(VEL**3)+0.508D0*(VEL**2)
     &    + 0.3755D0*(VEL) + 0.1213D0
        ELSEIF (VEL.LE.0.2D0) THEN
         INDV=0.2165D0
        ELSE
         INDV=0.9818D0
        ENDIF
!   SUMATORY OF INDEX
!    IND = SQRT( EXP(-2*(VEL - 1.0D0)**2) * EXP(-2*(H%R(JJ) - 1.0D0)**2))
        ECOSUI%ADR(II)%P%R(JJ) = SQRT(INDH*INDV)
        INDH = 0D0
        INDV = 0D0
!
       ENDDO
      ENDDO
!WRITE(*,*) "IIIIIIIND", VEL
!
!
! COMPUTING GRADIENT OF THE HABITAT SUITABILITY FIELD
!
      IELM = H%ELM
      CALL VECTOR(ECGRXN%ADR(1)%P,'=','GRADF          X',IELM,
     &              1D0,ECOSUI%ADR(1)%P,S,S,S,S,S,MESH,MSK,MASKEL)
!
      DO II = 1,NPOIN
        ECGRXN%ADR(1)%P%R(II) = ECGRXN%ADR(1)%P%R(II)*UNSV2D%R(II)
      END DO
! http://opentelemac.com/index.php/assistance/forum5/16-telemac-2d/8028-slope

      RETURN
      END
