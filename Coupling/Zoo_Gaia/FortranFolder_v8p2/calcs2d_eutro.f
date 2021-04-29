!                      ************************
                       MODULE DECLARATIONS_SUIT
!                      ************************
!
!***********************************************************************
! TELEMAC2D   V8P2
!***********************************************************************
!
!brief    DECLARATION OF PRINICIPAL ADDED VARIABLES
!
!history  Sergio Castiblanco
!         04/21
!+
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF_DEF
      USE FRICTION_DEF
!-----------------------------------------------------------------------
!
!       VARIABLES FOR COUPLING PURPOUSES
!
!-----------------------------------------------------------------------
!
!      NUMBER UF USER ADITIONAL VARS
!
      INTEGER :: NUSRVAR = 5
!
!     NEW MODIFIED VARIABLE (SOME INDEX)
!
      TYPE(BIEF_OBJ), TARGET :: IDX
!
!     BIEF_OBJ WITH THE COUECO NUMERATION (VECTOR)
!
!      TYPE(BIEF_OBJ), TARGET :: COUECO
!
!     BIEF_OBJ WITH THE FUNCTIONAL GROUPS RESULTS (BLOCK)
!
!      TYPE(BIEF_OBJ), TARGET :: ECOOUT
      TYPE(BIEF_OBJ), TARGET :: GRUF1
!
!     BIEF_OBJ WITH HABITAT SUITABILITY (BLOCK)
!
!      TYPE(BIEF_OBJ), TARGET :: ECOSUI
!
!
!     POINTERS TO FUNCTIONAL GROUPS RESULTS
!
      !TYPE(BIEF_OBJ), POINTER :: GRUF1, GRUF2, GRUF3, GRUF4, GRUF5
      !TYPE(BIEF_OBJ), POINTER :: GRUF6, GRUF7, GRUF8, GRUF9, GRUF10                 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SAVE
!
!     ============================================
!
      END MODULE DECLARATIONS_SUIT
!
!-----------------------------------------------------------------------
!
!        **********************************************
                 MODULE DECLARATIONS_ECOPLANKTON
!        **********************************************
!
!
!***********************************************************************
! WAQTEL   V8P2
!***********************************************************************
!
!brief    DECLARATION OF PRINCIPAL ECOPLANKTON ADDED VARIABLES TO WAQTEL
!
!history  Sergio Castiblanco
!+        04/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF_DEF
!
!       NOTE: THIS MODULE IS ORGANIZED IN 6 PARTS
!
!       (1) VECTORS (WILL BE DECLARED AS BIEF_OBJ STRUCTURES)
!       (2) INTEGERS
!       (3) LOGICAL VALUES
!       (4) REALS
!       (5) STRINGS
!       (6) ALIASES (IF NECESSARY)
!
!-----------------------------------------------------------------------
! (1) VECTORS (REAL AND INTEGER)
!-----------------------------------------------------------------------

!     FEEDING RATE OF THE ZOOPLANKTON FROM PHYTOPLANKTON, ZOOPLANKTON AN
!     D ORGANIC LOAD
      TYPE(BIEF_OBJ), TARGET :: ZOOGS1
      TYPE(BIEF_OBJ), TARGET :: ZOOGS2
      TYPE(BIEF_OBJ), TARGET :: ZOOGS3
!
!     WORKING ARRAYS FOR DEAL WITH ZOOPLANKTON
      TYPE(BIEF_OBJ), TARGET :: ZOOT1
      TYPE(BIEF_OBJ), TARGET :: ZOOT2
!
!-----------------------------------------------------------------------
!
!       2) INTEGERS
!
!-----------------------------------------------------------------------
! THIS IS THE IDENTIFIER OF ZOOPLANKTON INTO THE TRACERS BLOCK, TELEMAC
! ASSIGN IT A VALUE AUTOMATICALLY BY CALLING THE SUBROUTINE ADDTRACER
! INTO THE SUBROUTINE NAMETRAC_WAQTEL.F
      INTEGER IND_ZOO      
!
!-----------------------------------------------------------------------
!
!       3) LOGICAL
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!       4) REALS
!
!-----------------------------------------------------------------------
! PART ADDED BY SERGIO
!     ZOOPLANKTON CONSTANTS AND DUE TO ZOOPLANKTON CONSTANTS
!
!     SPECIFIC MORTALITY RATE OF THE ZOOPLANKTON (DAY^-1)
      DOUBLE PRECISION, PARAMETER :: MUZ = 0.1D0
!      DOUBLE PRECISION, PARAMETER :: MUZ = 0.0D0
!
!     SPECIFIC RATE OF ZOOPLANKTON METABOLIC EXCRETIONS (DAY^-1)
      DOUBLE PRECISION, PARAMETER :: GAMMAZ = 0.01D0
!      DOUBLE PRECISION, PARAMETER :: GAMMAZ = 0.0D0
!
!     SEMI-SATURATION CONSTANT BY FOOD IN THE ZOOPLANKTON GRWOTH PROCESS
!     (mgC m^-3)
      DOUBLE PRECISION, PARAMETER :: ZOOBK = 4150D0
!
!     MAXIMUM GRAZING RATE (DAY^-1)
      DOUBLE PRECISION, PARAMETER :: ZOOG = 0.75D0
!      DOUBLE PRECISION, PARAMETER :: ZOOG = 0D0
!
!     STOICHIOMETRIC TRANSFER COEFFICIENT FROM mgC TO mgP, FOR PO4
!     (mgC/mgP)
      DOUBLE PRECISION, PARAMETER :: BETAPC = 0.024D0
!      DOUBLE PRECISION, PARAMETER :: BETAPC = 0D0
!
!     TRANSFER COEFFICIENT FROM CUBIC METERS TO LITER, FOR VARIOUS
      DOUBLE PRECISION, PARAMETER :: BETAM3L = 0.001D0
!
!     STOICHIOMETRIC TRANSFER COEFFICIENT FROM mgC TO mgN, FOR NH4
!     (mgC/mgN)
      DOUBLE PRECISION, PARAMETER :: BETANC = 0.176D0
!      DOUBLE PRECISION, PARAMETER :: BETANC = 0D0
!
!     STOICHIOMETRIC TRANSFER COEFFICIENT FROM mgC TO O2
!     (mgO2/mgN)
      DOUBLE PRECISION, PARAMETER :: BETAO2C = 2.67D0
!      DOUBLE PRECISION, PARAMETER :: BETAO2C = 0D0
!
!     ASSIMILATION FEEDING OF ZOOPLANKTON FROM PHYTOPLANKTON, ZOOPLANKTO
!     N AND ORGANIC LOAD
      DOUBLE PRECISION, DIMENSION(3) :: ZOOOME
!
!     ZOOPLANKTON DIET FROM PHYTOPLANKTON, ZOOPLANKTON AND
!     ORGANIC LOAD
      DOUBLE PRECISION, DIMENSION(3) :: ZOORHOS
!
!-----------------------------------------------------------------------
!
!       5) STRINGS
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      SAVE
!
      END MODULE DECLARATIONS_ECOPLANKTON
!
!-----------------------------------------------------------------------
!
!                     ************************
                      SUBROUTINE CALCS2D_EUTRO
!                     ************************
!
     &  (NPOIN,WATTEMP,TN,TEXP,TIMP,RAYEFF,HPROP,T1,T2,T3,T4,
     &   T5,T6,T7,T8,T9,T10,T11,T12,DEBUG,UN,VN)
!
!***********************************************************************
! WAQTEL   V8P1
!***********************************************************************
!
!brief    COMPUTES SOURCE TERMS FOR EUTRO WAQ PROCESS
!
!history  R. ATA
!+        21/09/2014
!+        V7P0
!+       CREATION (VOID)
!history  R. ATA
!+        21/09/2015
!+        V7P1
!+       REAL IMPLEMENTATION
!
!history  R. ATA
!+        21/03/2016
!+        V7P2
!+       IMPROVEMENT- REMOVE LOCAL DECLARATIONS
!+       AND ALLOCATIONS
!
!history  S.E. BOURBAN (HRW)
!+        07/06/2017
!+        V7P3
!+        Indexing tracer (IND_*) to avoid conflicting naming convention
!+        between user defined tracers, water quality processes and
!+        ice processes. Introduction of the array RANK_*.
!
!history  S.E. BOURBAN (HRW)
!+        25/09/2017
!+        V7P3
!+        TEXP and TIMP are now additive to account for a variety of
!+        of sources / sinks on a given TRACER
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DEBUG          |-->| IF NE.0 THEN DEBUG MODE
!| HPROP          |-->| WATER DEPTH
!| NPOIN          |-->| TOTAL NUMBER OF MESH NODES
!| RAYEFF         |-->| EFFECT OF SUNSHINE ON ALGAE GROWTH
!| T1,..,T12      |<--| WORKING STRUCTURES
!| TN             |-->| TRACER STRUCUTRE
!| TEXP           |<--| EXPLICIT SOURCE TERMES
!| WATTEMP        |-->| WATER TEMPERATURE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!-----------------------------------------------------------------------
!***********************************************************************
      USE BIEF
      USE INTERFACE_PARALLEL
      USE DECLARATIONS_WAQTEL,ONLY:CMAX,CTOXIC,IK,K520,O2SATU,K2,ZSD,
     &  I0,KPE,KP,KN,CMORALG,FORMCS,TRESPIR,PROPHOC,DTP,PRONITC,K22,
     &  WLOR,K360,K320,PERNITS,WPOR,WNOR,FORMK2,O2PHOTO,K120,O2NITRI,
     &  DEMBEN,MEXTINC,SECTODAY,IND_T,
     &  IND_PHY,IND_PO4,IND_POR,IND_NO3,IND_NOR,IND_NH4,IND_OL,IND_O2
      USE INTERFACE_WAQTEL, EX_CALCS2D_EUTRO => CALCS2D_EUTRO
! .___________.____.____.______________________________________________.
! !    NOM    !TYPE!MODE!                   ROLE                       !
! !___________!____!____!______________________________________________!
! !  U        ! TR ! D  ! VITESSE DE L'EAU                             !
! !  WPOR     ! R  !    ! VITESSE DE SEDIMENTATION DU PHOSPHORE ORGANIQ!
! !  WNOR     ! R  !    ! VITESSE DE SEDIMENTATION DE L AZOTE ORGANIQUE!
! !  CMAX     ! R  !    ! TAUX DE CROISSANCE ALGALE MAXIMUM A 20°C     !
! !  ZSD      ! R  !    ! PROFONDEUR DE SECCHI                         !
! !  KPE      ! R  !    ! COEF D EXTINCTION DU RAY SANS PHYTO          !
! !  BETA     ! R  !    ! COEF DE TURBIDITE VEGETALE                   !
! !  IK       ! R  !    ! PARAMETRE DE CALAGE DE LA FORMULE DE SMITH   !
! !  KP       ! R  !    ! CONSTANTE DE DEMI-SATURATION EN PHOSPHATE    !
! !  KN       ! R  !    ! CONSTANTE DE DEMI-SATURATION EN AZOTE        !
! !  ALPHA    ! R  !    ! COEF 1 DE TOXICITE DE L EAU POUR LES ALGUES  !
! !  ALPHA2   ! R  !    ! COEF 2 DE TOXICITE DE L EAU POUR LES ALGUES  !
! !  RP       ! R  !    ! TAUX DE RESP. DE LA BIOMASSE ALGALE A 20°C   !
! !  PROPHOC  ! R  !    ! PROP DE PHOSPHORE DANS LES CELLULES DU PHYTO !
! !  DTP      ! R  !    ! POURCENT DE PHOSPH DIRECT ASSIM DS PHY MORT  !
! !  K320     ! R  !    ! TAUX DE TRANSFORMATION DU POR EN PO4         !
! !  PRONITC  ! R  !    ! PROP D AZOTE DANS LES CELLULES DU PHYTO      !
! !  PERNITS  ! R  !    ! POURCENT D AZOTE DIRECT ASSIM DS PHY MORT    !
! !  K360     ! R  !    ! TAUX DE TRANSFORMATION DU NOR EN NO3         !
! !  M1       ! R  !    ! COEF 1 DE MORTALITE ALGALE A 20°C            !
! !  M2       ! R  !    ! COEF 2 DE MORTALITE ALGALE A 20°C            !
! !  WLOR     ! R  !    ! VITESSE DE SEDIMENTATION DE LA CHARGE ORGANIQ!
! !  K120     ! R  !    ! CINETIQUE DE DEGRADATION DE LA CHARGE ORGANIQ!
! !  K520     ! R  !    ! CINETIQUE DE NITRIFICATION                   !
! !  F        ! R  !    ! QTTE D O2 PRODUITE PAR PHOTOSYNTHESE         !
! !  N        ! R  !    ! QTTE D O2 CONSOMMEE PAR NITRIFICATION        !
! !  BEN      ! R  !    ! DEMANDE BENTHIQUE A 20°C                     !
! !  K2       ! R  !    ! COEFFICIENT DE REAERATION                    !
! !  FORMK2   ! E  !    ! FORMULE DE CALCUL DE K2                      !
! !  CS       ! R  !    ! CONCENTRATION DE SATURATION EN OXYG DE L'EAU !
! !  FORMCS   ! E  !    ! FORMULE DE CALCUL DE CS                      !
! !  RSW      ! R  !    ! COEFFICIENT DE REAERATION AUX SEUILS         !
! !  FORMRS   ! E  !    ! FORMULE DE CALCUL DE R                       !
! !  ARS      ! R  !    ! COEFFICIENT A DES FORMULES DE CALCUL DE R    !
! !  BRS      ! R  !    ! COEFFICIENT B DES FORMULES DE CALCUL DE R    !
! !  NBSEUI   ! E  !    ! NOMBRE DE SEUILS                             !
! !  XSEUI    ! TR !    ! ABSCISSES DES SEUILS                         !
! !  DZS      ! TR !    ! DELTA Z AUX SEUILS                           !
! !           !    !    !                                              !
! !  IF1      ! TR ! D  ! INDIC DE LECTURE DU FICHIER DES PARAMETRES   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDED PART BY SERGIO
! !  MUZ      ! R  !    ! SPECIFIC MORTALITY RATE OF THE ZOOPLAKNTON   !
! !  GAMMAZ   ! R  !    ! SPECIFIC RATE OF ZOOPLANKTON METABOLIC EXCR. !
! !  ZOOBK    ! R  !    ! SEMI-SATURATION CONSTANT BY FOOD IN ZOO GROWT!
! !  ZOOG     ! R  !    ! MAXIMUM ZOOPLANKTON GRAZING RATE             !
! !  BETAPC   ! R  !    ! TRANSFER COEFF. FROM mgC TO mgP              !
! !  BETAM3L  ! R  !    ! TRANSFER COEFF. FROM m^3 TO LITER            !
! !  BETANC   ! R  !    ! TRANSFER COEFF. FROM mgC TO mgN, FOR NH4     !
! !  BETAO2C  ! R  !    ! TRANSFER COEFFICIENT FROM mgC TO mgO2        !
! !  ZOOGS1,2.! R  !    ! FEEDING RATES OF ZOO FROM PHYTO, ZOO AND ORG.!
! !  ZOOPES   ! R  !    ! BIOMASS RELATIONS BETWEEN PHYTO, ZOO AND ORG.!
! !  ZOOOME   ! R  !    ! ASSIMILATION COEFF. OF ZOO FEEDING PROCESS   !
! !  ZOORHOS  ! R  !    ! ZOOPLANKTON DIET COMPOSITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !___________!____!____!______________________________________________!
!  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
!               (ENTREE)              (SORTIE)       (ENTREE/SORTIE)
!-----------------------------------------------------------------------
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART ADDED BY SERGIO
! USING DECLARATIONS_ECOPLANKTON
      USE DECLARATIONS_ECOPLANKTON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER          , INTENT(IN   ) :: NPOIN,DEBUG
      DOUBLE PRECISION , INTENT(IN   ) :: WATTEMP
      TYPE(BIEF_OBJ)   , INTENT(IN   ) :: TN,HPROP,UN,VN
      TYPE(BIEF_OBJ)   , INTENT(INOUT) :: TIMP,TEXP,RAYEFF
      TYPE(BIEF_OBJ)   , INTENT(INOUT) :: T1,T2,T3,T4,T5,T6,T7,T8,T9,T10
      TYPE(BIEF_OBJ)   , INTENT(INOUT) :: T11,T12
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     LOCAL VARIABLES
!
      INTEGER                     :: I
      DOUBLE PRECISION, PARAMETER :: UNSURVINGT=0.05D0
      DOUBLE PRECISION, PARAMETER :: EPS=1.D-6
      DOUBLE PRECISION, PARAMETER :: CORR1=1.065D0
!      DOUBLE PRECISION, PARAMETER :: CORR2=1.0241D0
      DOUBLE PRECISION, PARAMETER :: CORR2=1.025D0
      DOUBLE PRECISION            :: G1,G2,G3,CC,POWER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDED PART BY SERGIO
! DUMMY VALUE TO MAKE THINGS CLEAR
      DOUBLE PRECISION :: DUMMY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 0'
!
!     INITIALISATION
!
!     CS IS STORED IN T2
      CALL OS( 'X=0     ',X=T2)
      CALL OS( 'X=0     ',X=T3)
!     G2 IS STORED IN T6
      CALL OS( 'X=0     ',X=T6)
!     CP IS STORED IN T7
      CALL OS( 'X=0     ',X=T7)
!     DP IS STORED IN T4
      CALL OS( 'X=0     ',X=T4)
!     LNUT IS STORED IN T5 UNTIL ALGAE_GROWTH,
!     THEN RESPIR*G1 IS STORED IN T5 FOR THE 8TH EQUATION
      CALL OS( 'X=0     ',X=T5)
!     RN IS STORED IN T8
      CALL OS( 'X=0     ',X=T8)
!
!     G2 IS STORED IN T6,WE TAKE INTO ACCOUNT VARIABLE TEMPERATURE
!
!      G2 = WATTEMP/20.D0
!      IF( IND_T.GT.0 ) THEN
!        CALL OS('X=CY    ',X=T6,Y=TN%ADR(IND_T)%P,C=UNSURVINGT)
!      ELSE
!        CALL OS('X=C     ',X=T6 ,C=G2)
!      ENDIF
!
!     COMPUTE G3,BENCOR
!
      POWER = WATTEMP-20.D0
      G2    = 1.050D0**POWER
      G3    = 1.047D0**POWER
      DO I=1,NPOIN
        IF( IND_T.GT.0 ) THEN
          POWER = TN%ADR(IND_T)%P%R(I)-20.D0
          G2    = 1.050D0**POWER
          G3    = 1.047D0**POWER
        ENDIF
!       CORR2T AND BENCOR STORED HERE IN T9,T10
        T9%R(I) = CORR2**POWER
        T10%R(I)= DEMBEN*(CORR1**POWER)
!       G2 IS STORED IN T6
        T6%R(I)=G2
!       G3 IS STORED IN T11
        T11%R(I)=G3
      ENDDO
!
!     COMPUTE CS (O2SATU, STORED IN T2)
!
      IF( IND_T.EQ.0 ) THEN
        CALL SATUR_O2(O2SATU,FORMCS,WATTEMP,EPS)
        CALL OS('X=C     ',X=T2,C=O2SATU       )
      ELSE
        DO I=1,NPOIN
          CALL SATUR_O2(T2%R(I),FORMCS,TN%ADR(IND_T)%P%R(I),EPS)
        ENDDO
      ENDIF
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 1'
!
!     RAYEFF WITH SMITH FORMULA
!
      CALL RAY_EFFECT(ZSD,TN%ADR(IND_PHY)%P,NPOIN,MEXTINC,I0,IK,KPE,
     &                RAYEFF,HPROP,T3,T4)
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 2'
!
!     COMPUTE LNUT: EFFECTS OF PHOSPHORIOUS AND NITROGENIOUS
!           NUTRIMENTS ON ALGAE GROWTH ==>STORED IN T5
!
!      CALL NUTEFF(T5%R,TN,NPOIN,IND_PO4,IND_NO3,KP,KN)
!
!     NUTEFF DOEST NOT TO NH4 INTO ACCOUNT
      DO I=1,NPOIN
        T5%R(I)= MIN(TN%ADR(IND_PO4)%P%R(I)/(KP+TN%ADR(IND_PO4)%P%R(I)),
     &                (TN%ADR(IND_NO3)%P%R(I)+TN%ADR(IND_NH4)%P%R(I))
     &               /(KN+TN%ADR(IND_NO3)%P%R(I)
     &                   +TN%ADR(IND_NH4)%P%R(I)))
      ENDDO
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 3'
!
!     RATE OF ALGAE GROWTH: CP (STORED IN T7)
!
      CALL ALGAE_GROWTH(T7%R,CMAX,RAYEFF%R,T6,T5%R,CTOXIC(1),NPOIN)
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 4'
!
!     RATE OF ALGAE DISAPPEARANCE DP (STORED IN T4) AND MP (STOPCKED IN T12)
!
      CALL ALGAE_DEATH(T4%R,T12%R,CMORALG,TN%ADR(IND_PHY)%P%R,TRESPIR,
     &                  T6,CTOXIC(2),NPOIN)
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 5'
!
!     COMPUTE K2
!
      CALL REAER(FORMK2,K2,K22,NPOIN,1,UN,VN,HPROP,EPS)
!
!     COMPUTE RS (RSW:DONE IN DIFSOU)
!
!
!     COMPUTE RN: PROPORTION OF NITROGEN ASSIMILATED AS NH4(STORED IN T8)
!
      CALL OV( 'X=Y+Z   ' ,T1%R,TN%ADR(IND_NH4)%P%R,TN%ADR(IND_NO3)%P%R,
     &          0.D0,NPOIN)
      CALL OVD('X=Y/Z   ' ,T8%R,TN%ADR(IND_NH4)%P%R,T1%R,0.D0,
     &          NPOIN ,2,0.D0,EPS )
!
!     RESPIR*G1 IS STORED IN T5 FOR THE 8TH EQUATION
      G1 = WATTEMP/20.D0
      IF( IND_T.GT.0 ) THEN
        CALL OS('X=CY    ',X=T5,Y=TN%ADR(IND_T)%P,C=UNSURVINGT*TRESPIR)
      ELSE
        CALL OS('X=C     ',X=T5 ,C=G1*TRESPIR)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDED PART BY SERGIO
!
! ZOOPLANKTON CALCULATIONS ARE COMPUTED HERE
      CALL ZOO_CALCS(ZOOT1,ZOOT2,ZOOGS1,ZOOGS2,ZOOGS3,ZOOBK,ZOORHOS,
     & ZOOOME, TN%ADR(IND_PHY)%P,TN%ADR(IND_OL)%P,TN%ADR(IND_ZOO)%P,
     & NPOIN, GAMMAZ,MUZ,BETAM3L,ZOOG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!
!     LET'S NOW COMPUTE SOURCE TERMS
!     -------------------------------
!
!     FIRST TRACER [PHY] (IND_PHY)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED PART BY SERGIO
! ORIGINAL
!      CALL OS( 'X=Y-Z   ' ,X=T1                 ,Y=T7,Z=T4)
!!      CALL OS( 'X=YZ    ' ,X=TEXP%ADR(IND_PHY)%P,Y=T1,
!!     &                     Z=TN%ADR(IND_PHY)%P)
!      CALL OS( 'X=X+CYZ ' ,X=TEXP%ADR(IND_PHY)%P,Y=T1,
!     &                     Z=TN%ADR(IND_PHY)%P,C=SECTODAY )
!
! MODIFIED
! HERE AND UNTIL THE END OF THIS SUBROUTINE, THE CHANGES DUE TO ZOOPLAN.
! INCLUSION ARE MADE
!  CP - DP
      CALL OS( 'X=Y-Z   ' ,X=T1                 ,Y=T7,Z=T4)
!  (CP - DP)*PHY
      CALL OS( 'X=XY    ' ,X=T1,Y=TEXP%ADR(IND_PHY)%P)
!  (CP - DP)*PHY - Gs1*ZOO
      CALL OS( 'X=X-YZ  ' ,X=T1,Y=ZOOGS1,Z=TN%ADR(IND_ZOO)%P        )
!  ((CP - DP)*PHY - Gs1*ZOO)*SECTODAY
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_PHY)%P,Y=T1,C=SECTODAY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 6'
!
!     SECOND TRACER [PO4] (IND_PO4)
!
!  dtp*DP
      CALL OS( 'X=CY    ' ,X=T1,Y=T4                    ,C=DTP      )
!  dtp*DP - CP
      CALL OS( 'X=X-Y   ' ,X=T1,Y=T7                                )
!  fp(dtp*DP - CP)*PHY
      CALL OS( 'X=CXY   ' ,X=T1,Y=TN%ADR(IND_PHY)%P     ,C=PROPHOC  )
!  k320*g2*POR
      CALL OS( 'X=CYZ   ' ,X=T3,Y=TN%ADR(IND_POR)%P,Z=T6,C=K320     )
!      CALL OS( 'X=Y+Z   ' ,X=TEXP%ADR(IND_PO4)%P,Y=T1,Z=T3          )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED PART BY SERGIO
! ORIGINAL
!      CALL OS( 'X=X+Y   ' ,X=T1,Y=T3                                )
!      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_PO4)%P,Y=T1   ,C=SECTODAY )
! MODIFIED
      DUMMY = GAMMAZ*BETAPC*BETAM3L
      CALL OS( 'X=X+CY  ' ,X=T3,Y=TN%ADR(IND_ZOO)%P     ,C=DUMMY    )
      CALL OS( 'X=X+Y   ' ,X=T1,Y=T3                                )
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_PO4)%P,Y=T1   ,C=SECTODAY )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 7'
!
!     THIRD TRACER [POR] (IND_POR)
!
      G2=PROPHOC*(1.D0-DTP)
      CALL OS( 'X=CYZ   ' ,X=T1,Y=T4,Z=TN%ADR(IND_PHY)%P,C=G2       )
      CALL OS( 'X=X-Y   ' ,X=T1,Y=T3                                )
!      CALL OVD('X=C/Y   ' ,T3%R,HPROP%R,TN%ADR(IND_POR)%P%R,WPOR,
!     &          NPOIN ,2,0.D0,EPS                                   )
      CALL OVD('X=CY/Z  ' ,T3%R,TN%ADR(IND_POR)%P%R,HPROP%R,WPOR,
     &         NPOIN ,2,0.D0,EPS )
!      CALL OS( 'X=Y-Z   ' ,X=TEXP%ADR(IND_POR)%P,Y=T1,Z=T3          )
      CALL OS( 'X=X-Y   ' ,X=T1,Y=T3                                )
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_POR)%P,Y=T1   ,C=SECTODAY )
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 8'
!
!     FOURTH TRACER [NO3] (IND_NO3)
!
      CALL OS( 'X=Y+C   ' ,X=T1,Y=T8                      ,C=-1.D0  )
      CALL OS( 'X=CXYZ  ' ,X=T1,Y=T7,Z=TN%ADR(IND_PHY)%P  ,C=PRONITC)
      CALL OS( 'X=CYZ   ' ,X=T3,Y=TN%ADR(IND_NH4)%P  ,Z=T6,C=K520)
!      CALL OS( 'X=Y-Z   ' ,X=TEXP%ADR(IND_NO3)%P,Y=T3,Z=T1          )
      CALL OS( 'X=X+Y   ' ,X=T3,Y=T1                                )
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_NO3)%P,Y=T3   ,C=SECTODAY )
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 9'
!
!     FIFTH TRACER [NOR] (IND_NOR)
!
      G2=PRONITC*(1.D0-PERNITS)
      CALL OS( 'X=CYZ   ' ,X=T1,Y=T4,Z=TN%ADR(IND_PHY)%P,C=G2       )
      CALL OS( 'X=CYZ   ' ,X=T3,Y=TN%ADR(IND_NOR)%P  ,Z=T6,C=K360)
      CALL OS( 'X=X-Y   ' ,X=T1,Y=T3                                )
      CALL OVD('X=CY/Z  ' ,T3%R,TN%ADR(IND_NOR)%P%R,HPROP%R,WNOR,
     &         NPOIN ,2,0.D0,EPS )
!      CALL OS( 'X=Y-Z   ' ,X=TEXP%ADR(IND_NOR)%P,Y=T1,Z=T3          )
      CALL OS( 'X=X-Y   ' ,X=T1,Y=T3                                )
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_NOR)%P,Y=T1   ,C=SECTODAY )
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 10'
!
!     SIXTH TRACER [NH4] : AMMONIACAL LOAD (IND_NH4)
!
!     IMPLICIT PART
!      CALL OS( 'X=CYZ   ' ,X=TIMP%ADR(IND_NH4)%P,Y=T6,Z=HPROP,C=-K520)
      CALL OS( 'X=X+CYZ ' ,X=TIMP%ADR(IND_NH4)%P,
     &                     Y=T6,Z=HPROP,   C=-K520*SECTODAY)
!     EXPLICIT PART
!  dtn*DP
      CALL OS( 'X=CY    ' ,X=T1,Y=T4                      ,C=PERNITS )
!  Rn*CP
      CALL OS( 'X=YZ    ' ,X=T3,Y=T7,Z=T8                            )
!  fn*(dtn*DP - Rn*CP)
      CALL OS( 'X=C(Y-Z)' ,X=T3,Y=T1,Z=T3                 ,C=PRONITC )
!  fn*(dtn*DP - Rn*CP)*PHY
      CALL OS( 'X=XY    ' ,X=T3,Y=TN%ADR(IND_PHY)%P                  )
!  k620*g2*NOR      TYPO: k620 Theory is k360 Code
      CALL OS( 'X=CYZ   ' ,X=T1,Y=TN%ADR(IND_NOR)%P  ,Z=T6,C=K360    )
!      CALL OS( 'X=Y+Z   ' ,X=TEXP%ADR(IND_NH4)%P,Y=T1,Z=T3           )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED PART BY SERGIO
! ORIGINAL
!      CALL OS( 'X=X+Y   ' ,X=T1,Y=T3                                )
!      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_NH4)%P,Y=T1   ,C=SECTODAY )
! MODIFIED
      DUMMY = GAMMAZ*BETANC*BETAM3L
      CALL OS( 'X=X+CY  ' ,X=T1,Y=TN%ADR(IND_ZOO)%P     ,C=DUMMY    )
      CALL OS( 'X=X+Y   ' ,X=T1,Y=T3                                )
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_NH4)%P,Y=T1   ,C=SECTODAY )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 11'
!
!     SEVENTH TRACER [L]: ORGANIC LOAD (IND_OL)
!
!     IMPLICIT PART
!      CALL OS( 'X=CYZ   ' ,X=TIMP%ADR(IND_OL)%P,Y=T11,Z=HPROP,C=-K120)
      CALL OS( 'X=X+CYZ ' ,X=TIMP%ADR(IND_OL)%P,
     &                     Y=T11,Z=HPROP,C=-K120*SECTODAY )
!     EXPLICIT PART
!  f*MP*PHY
      CALL OS( 'X=CYZ   ' ,X=T1,Y=T12,Z=TN%ADR(IND_PHY)%P,C=O2PHOTO  )
!  FLOR/H
      CALL OVD('X=CY/Z  ' ,T3%R,TN%ADR(IND_OL)%P%R,HPROP%R,WLOR,
     &          NPOIN ,2,0.D0,EPS                                    )
!      CALL OS( 'X=Y-Z   ' ,X=TEXP%ADR(IND_OL)%P,Y=T1,Z=T3            )
      CALL OS( 'X=X-Y   ' ,X=T1,Y=T3                                 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED PART BY SERGIO
! ORIGINAL
!      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_OL)%P,Y=T1   ,C=SECTODAY   )
! MODIFIED
      CALL OS( 'X=X+CYZ ' ,X=T1,Y=ZOOT2,Z=TN%ADR(IND_ZOO)%P,C=BETAM3L)
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_OL)%P,Y=T1   ,C=SECTODAY   )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 12'
!
!     EIGTH TRACER: DISSOLVED O2 (IND_O2)
!
!      CALL OS( 'X=Y+C   ' ,X=T1                 ,Y=T7,C=-TRESPIR     )
!  (CP-RP*g1)
      CALL OS( 'X=Y-Z   ' ,X=T1                 ,Y=T7,Z=T5           )
!      CALL OS( 'X=CYZ   ' ,X=TEXP%ADR(IND_O2)%P,Y=T1,
!  f*(CP-RP*g1)*PHY
      CALL OS( 'X=CYZ   ' ,X=T4,Y=T1,
     &                     Z=TN%ADR(IND_PHY)%P,C=O2PHOTO)
!  n*k520*g2*NH4
!     -nK520g2[NH4]
      CC=O2NITRI*K520
      CALL OS( 'X=CYZ   ' ,X=T1,Y=T6,Z=TN%ADR(IND_NH4)%P,C=CC        )
!      CALL OS( 'X=X-Y   ' ,X=TEXP%ADR(IND_O2)%P,Y=T1                 )
!  f*(CP-RP*g1)*PHY - n*k520*g2*NH4
      CALL OS( 'X=X-Y   ' ,X=T4,Y=T1                                 )
!  k120*g3*L
!     K120g3[L]
      CALL OS( 'X=CYZ   ' ,X=T1,Y=T11,Z=TN%ADR(IND_OL)%P,C=K120)
!      CALL OS( 'X=X-Y   ' ,X=TEXP%ADR(IND_O2)%P,Y=T1                 )
!  f*(CP-RP*g1)*PHY - n*k520*g2*NH4 - k120*g3*L
      CALL OS( 'X=X-Y   ' ,X=T4,Y=T1                                 )
!     K2g4(Cs-[O2])
!  Cs - O2
      CALL OS( 'X=Y-Z   ' ,X=T1,Y=T2,Z=TN%ADR(IND_O2)%P              )
!  k2*g4*(Cs - O2)
      CALL OS( 'X=CXYZ  ' ,X=T1,Y=T9,Z=K2,C=1.D0                     )
!      CALL OS( 'X=X+Y   ' ,X=TEXP%ADR(IND_O2)%P,Y=T1                 )
!  f*(CP-RP*g1)*PHY - n*k520*g2*NH4 - k120*g3*L + k2*g4*(Cs - O2)
      CALL OS( 'X=X+Y   ' ,X=T4,Y=T1                                )
!     -BEN/h
      CALL OVD('X=Y/Z   ' ,T3%R,T10%R,HPROP%R,0.D0,
     &          NPOIN ,2,0.D0,EPS )
!      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_O2)%P,Y=T3,C=-DEMBEN       )
!      CALL OS( 'X=X+CY  ' ,X=T4,Y=T3                ,C=-DEMBEN       )
!  f*(CP-RP*g1)*PHY - n*k520*g2*NH4 - k120*g3*L + k2*g4*(Cs - O2) - BEN/h
      CALL OS( 'X=X-Y   ' ,X=T4,Y=T3                                 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED PART BY SERGIO
! ORIGINAL
!      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_O2)%P,Y=T4,C=SECTODAY      )
! MODIFIED
      DUMMY = -1D0*GAMMAZ*BETAO2C*BETAM3L
      CALL OS( 'X=X+CY  ' ,X=T4,Y=TN%ADR(IND_ZOO)%P,C=DUMMY)
      CALL OS( 'X=X+CY  ' ,X=TEXP%ADR(IND_O2)%P,Y=T4,C=SECTODAY      )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDED PART BY SERGIO
!
!     NINETH TRACER: ZOOPLANKTON (IND_ZOO)
!
      CALL OS( 'X=X+CYZ ' ,X=TEXP%ADR(IND_ZOO)%P,Y=ZOOT1,
     &                     Z=TN%ADR(IND_ZOO)%P  ,C=SECTODAY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(DEBUG.GT.0)WRITE(LU,*)'IN EUTRO, STEP 14'
!
!-----------------------------------------------------------------------
!
      RETURN
      END
