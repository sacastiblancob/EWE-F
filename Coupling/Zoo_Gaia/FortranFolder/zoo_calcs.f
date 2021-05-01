!                    **********************
                     SUBROUTINE ZOO_CALCS
!                    **********************
!
     & (ZOOT1,ZOOT2,ZOOGS1,ZOOGS2,ZOOGS3,ZOOBK,ZOORHOS,ZOOOME,
     &  PHYB,OLB,ZOOB,NPOIN,GAMMAZ,MUZ,ZOOG)
!
!***********************************************************************
! TELEMAC2D   V8P1R0
!***********************************************************************
!
!brief    COMPUTES ZOOPLANKTON DIETS RATES
!
!history  CASTIBLANCO-BALLESTEROS S. A.
!+        30/11/2020
!+        V8P1
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ZOOT1          |<--| TOTAL RATE COEFFICIENT FOR ZOOPLANCTON
!| ZOOT2          |<--| TOTAL RATE COEFFICIENT FOR ORG.L. EQUATION
!| ZOOGS1,2,3     |<->| PREDATION OF ZOOPLANKTON OVER PHYTO, ZOO AND OL
!| ZOOBK          |-->| SEMI-SATURATION CONSTANT TO COMPUTE ZOOGS
!| ZOORHOS        |-->| DIET COMPOSITION OF ZOOPLANKTON
!| ZOOOME         |-->| ASSIMILATION COEFFICIENTS OF ZOOPLANKTON
!| NPOIN          |-->| TOTAL NUMBER OF MESH NODES
!| PHYB           |-->| PHYTOPLANKTON
!| OLB            |-->| ORGANIC LOAD
!| ZOOB           |-->| ZOOPLANKTON
!| GAMMAZ         |-->| SPECIFIC RATE OF ZOOPLANKTON METABOLIC EXCR.
!| MUZ            |-->| SPECIFIC MORTALITY RATE OF THE ZOOPLAKNTON
!| ZOOG           |-->| MAXIMUM ZOOPLANKTON GRAZING RATE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
!
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,          INTENT(IN) :: NPOIN
      DOUBLE PRECISION, INTENT(IN) :: ZOOBK,GAMMAZ,MUZ,ZOOG
      DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: ZOORHOS,ZOOOME
      TYPE(BIEF_OBJ),   INTENT(IN) :: PHYB,OLB,ZOOB
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ZOOGS1,ZOOGS2,ZOOGS3
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ZOOT1,ZOOT2
! LOCAL VARIABLES
      INTEGER                        :: KK
      DOUBLE PRECISION               :: RHOB, SUMPB
      DOUBLE PRECISION, DIMENSION(3) :: ZOOPES
      DOUBLE PRECISION               :: EPSIL = 1E-10
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     INITIALISATION
!
      CALL OS( 'X=0     ' ,X=ZOOT1)
      CALL OS( 'X=0     ' ,X=ZOOT2)
      CALL OS( 'X=0     ' ,X=ZOOGS1)
      CALL OS( 'X=0     ' ,X=ZOOGS2)
      CALL OS( 'X=0     ' ,X=ZOOGS3)
      ZOOPES = 0D0
!
! COMPUTING FEEDING RATES WITH THE MODEL BASED ON ZOOPLANKTON DIET
!
      DO KK=1,NPOIN
        IF (ZOOB%R(KK) .LE. EPSIL) THEN
          ZOOGS1%R(KK) = 0D0
          ZOOGS2%R(KK) = 0D0
          ZOOGS2%R(KK) = 0D0
        ELSE
          ! COMPUTING SUMMATION OF BIOMASi TIMES RHOi
          RHOB = ZOORHOS(1)*PHYB%R(KK) + ZOORHOS(2)*ZOOB%R(KK) +
     &      ZOORHOS(3)*OLB%R(KK)
          ! COMPUTING BIOMASS RELATIONS Pj OR Pk
          ZOOPES(1) = ZOORHOS(1)*PHYB%R(KK)/RHOB
          ZOOPES(2) = ZOORHOS(2)*ZOOB%R(KK)/RHOB
          ZOOPES(3) = ZOORHOS(3)*OLB%R(KK)/RHOB
          ! COMPUTING SUMMATION OF Pj's OR Pk's
          SUMPB = ZOOPES(1) + ZOOPES(2) + ZOOPES(3)
          ! COMPUTING FEEDING RATES
          ZOOGS1%R(KK) = ZOOG*ZOORHOS(1)*PHYB%R(KK)/(ZOOBK + SUMPB)
          ZOOGS2%R(KK) = ZOOG*ZOORHOS(2)*ZOOB%R(KK)/(ZOOBK + SUMPB)
          ZOOGS3%R(KK) = ZOOG*ZOORHOS(3)*OLB%R(KK)
          ZOOGS3%R(KK) = ZOOGS3%R(KK)/(ZOOBK + SUMPB)
        ENDIF
      ENDDO
!
! COMPUTING FINAL RATE COEFFICIENT FOR ZOOPLANKTON EQUATION, STORED IN
!   ZOOT1
!
      DO KK=1,NPOIN
        IF (ZOOB%R(KK) .LE. EPSIL) THEN
          ZOOT1%R(KK) = 0D0
        ELSE
          ZOOT1%R(KK) = ZOOOME(1)*ZOOGS1%R(KK) + 
     &     (ZOOOME(2) - 1D0)*ZOOGS2%R(KK) + ZOOOME(3)*ZOOGS3%R(KK) -
     &     GAMMAZ - MUZ
        ENDIF
      ENDDO
!
! COMPUTING FINAL RATE COEFFICIENT FOR ORGANIC LOAD DUE TO ZOOPLANKTON
!   STORED IN ZOOT2
!
      DO KK=1,NPOIN
        IF (ZOOB%R(KK) .LE. EPSIL) THEN
          ZOOT2%R(KK) = 0D0
        ELSE
          ZOOT2%R(KK) = ((1D0 - ZOOOME(1))*ZOOGS1%R(KK) +
     &     (1D0 - ZOOOME(2))*ZOOGS2%R(KK) - ZOOOME(3)*ZOOGS3%R(KK) + 
     &      MUZ)
        ENDIF
      ENDDO
!   
!
      RETURN
      END



