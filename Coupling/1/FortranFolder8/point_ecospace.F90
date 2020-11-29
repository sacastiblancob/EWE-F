!                    **************************
                     SUBROUTINE POINT_ECOSPACE &
!                    **************************
  (nvars)
!***********************************************************************
! TELEMAC2D WITH ECOSPACE
!***********************************************************************
!
! Brief. Allocates all the BIEF_OBJ for coupling with ecospace
!     + created by Sergio Castiblanco - 26/06/2020
!
!***********************************************************************

      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE DECLARATIONS_SPECIAL
!      USE DECLARATIONS_COUPECOSPACE

!      USE statevartypesecosim, only: nvars

      IMPLICIT NONE

      !INTENT VARIABLES
      INTEGER, INTENT(IN)  :: nvars

      !IN SUBROUTINE VARIABLES
      INTEGER    :: II
      
      !ALLOCATES COUECO VECTOR
      CALL BIEF_ALLVEC(1,COUECO,'COUECO',IELMT,1,2,MESH)

      !ALLOCATES FUNCTIONAL GROUP RESULTS
      CALL ALLBLO(ECOOUT,'ECOOUT')
      CALL BIEF_ALLVEC_IN_BLOCK(ECOOUT,nvars,1,'GRUF  ',IELMT,1,2,MESH)

      CALL BIEF_ALLVEC(1,GRUF1,'GRUF1 ',IELMT,1,2,MESH)
!      GRUF1 => ECOOUT%ADR(1)%P

      !ALLOCATES HABITAT SUITABILITY BLOCK
      CALL ALLBLO(ECOSUI, 'ECOSUI')
      CALL BIEF_ALLVEC_IN_BLOCK(ECOSUI,nvars,1,'ESUI  ',IELMT,1,2,MESH)

!!!!!!!!!! INITILIASING BLOCKS IN ZERO !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      DO II = 1, nvars
          CALL OS('X=C     ',X=ECOOUT%ADR(II)%P,C=0.D0)
          CALL OS('X=C     ',X=ECOSUI%ADR(II)%P,C=0.D0)
      ENDDO
!      WRITE(*,*) "BBBBBB", ECOOUT%ADR(1)%P%R(1) 

      RETURN
      END