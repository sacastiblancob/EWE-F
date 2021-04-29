!        **********************************************
                   MODULE DECLARATIONS_WAQTEL
!        **********************************************
!
!
!***********************************************************************
! WAQTEL   V7P3
!***********************************************************************
!
!brief    DECLARATION OF PRINCIPAL WAQTEL VARIABLES
!
!history  R. ATA (EDF-LNHE)
!+
!+
!history  M.JODEAU (EDF-LNHE)
!+        08/2016
!+        V7P3
!+        AED2 coupling
!
!history  S.E. BOURBAN (HRW)
!+        07/06/2017
!+        V7P3
!+        Indexing tracer (IND_*) to avoid conflicting naming convention
!+        between user defined tracers, water quality processes and
!+        ice processes.
!
!history  S.E. BOURBAN (HRW)
!+        21/09/2017
!+        V7P3
!+        WAQPROCESS is now a prime number, so that multiple processes
!+        can be called by multiplication of the prime numbers.
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
!
!
!     COEFFICIENT OF REAERATION
!
      TYPE(BIEF_OBJ), TARGET :: K2
!
!     EFFECT OF SUNSHINE ON ALGAE GROWTH (BIOMASS PROCESS)
!
      TYPE(BIEF_OBJ), TARGET :: RAYEFF
!-----------------------------------------------------------------------
!
!       2) INTEGERS
!
!-----------------------------------------------------------------------
!
!     MAXIMUM NUMBER OF VALUES OF KEY-WORDS OF ONE KIND (INTEGER, LOGICAL, ETC.)
!
      INTEGER, PARAMETER :: MAXKEY = 300
!     GRAPHIC PRINTOUT PERIOD WAQ
!
      INTEGER LEOPRD
!
!     LISTING PRINTOUT PERIOD WAQ
!
      INTEGER LISPRD
!
!     MAXIMUM NUMBER OF OUTPUT VARIABLES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! PART MODIFIED BY SERGIO
!ORIGINAL      INTEGER, PARAMETER :: MAXWQVAR = 50
! MAXIMUM WAQTEL VARIABLES
      INTEGER, PARAMETER :: MAXWQVAR = 51
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------------------
!
!     WAQ PROCESS, A MULTPLICATIVE COMBINATION OF PRIME NUMBERS:
!     -  2: O2,
!     -  3: BIOMASS,
!     -  5: EUTRO,
!     -  7: MICROPOL,
!     - 11: THERMIC
!     - 13: AED2 LIBRARY
!     - 17: TRACER DEGRADATION
!     - 19: GHOST PROCESS IN WAITING FOR THE MERGE WITH ICE PROCESS
!
      INTEGER WAQPROCESS
!
!     TOTAL NUMBER OF TRACERS RELATED TO WAQPROCESS
!       AND THEIR ASSOCIATED INDICES
!
      INTEGER WAQTR
      INTEGER RANKTR(MAXWQVAR)
!
!     TOTAL NUMBER OF TRACERS RELATED TO INDIVIDUAL PROCESSES
!       AND THEIR ASSOCIATED INDICES
!
      INTEGER, PARAMETER :: NWAQ_OXYGN = 3
      INTEGER RANK_OXYGN(NWAQ_OXYGN)
!
      INTEGER, PARAMETER :: NWAQ_BIOMA = 5
      INTEGER RANK_BIOMA(NWAQ_BIOMA)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! PART MODIFIED BY SERGIO
!ORIGINAL
!      INTEGER, PARAMETER :: NWAQ_EUTRO = 8
!MODIFIED
! NUMBER OF TRACERS IN WAQTEL
      INTEGER, PARAMETER :: NWAQ_EUTRO = 9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER RANK_EUTRO(NWAQ_EUTRO)
!
      INTEGER, PARAMETER :: NWAQ_MICRO = 5
      INTEGER RANK_MICRO(NWAQ_MICRO)
!
      INTEGER, PARAMETER :: NWAQ_THERM = 1
      INTEGER RANK_THERM(NWAQ_THERM)
!
      INTEGER :: NWAQ_AED2
      INTEGER RANK_AED2(MAXWQVAR)
!
      INTEGER :: NWAQ_DEGRA
      INTEGER RANK_DEGRA(MAXWQVAR)
!
!     TO KNOW IF A VARIABLE WILL BE EXITED ON FILE, ON LISTING
!
      LOGICAL SORLEO(MAXWQVAR),SORIMP(MAXWQVAR)
!
!     WAQ RESULT FILE NUMBER
!
      INTEGER WAQRES
!
!     WAQ BOUNDARY CONDITION FILE NUMBER
!
      INTEGER WAQCLI
!
!     WAQ GEOMETRY FILE NUMBER
!
      INTEGER WAQGEO
!
!     WAQ HYDRODYNAMICS FILE NUMBER
!
      INTEGER WAQHYD
!
!     WAQ STEERING FILE NUMBER
!
      INTEGER WAQCAS
!
!     WAQ REFERENCE FILE NUMBER
!
      INTEGER WAQREF
!
!     DEBUGGER
!
      INTEGER DEBUG
!
!     FORMULA FOR COMPUTING K2
!
      INTEGER FORMK2
!
!     FORMULA FOR COMPUTING RS
!
      INTEGER FORMRS
!
!     FORMULA FOR COMPUTING CS
!
      INTEGER FORMCS
!
!     COEFFICIENTS OF AERATION FORMULA
!
      INTEGER CFORMAERA(2)
!
!     ATMOSPHERE-WATER EXCHANGE MODEL
!
      INTEGER ATMOSEXCH
!
!     BRIGHTNESS OF THE SKY
!
      INTEGER ISKYTYPE
!
!     METHOD OF COMPUTATION OF EXTINCTION OF SUN RAY
!
      INTEGER MEXTINC
!
!     FORMULA OF ATMOSPHERIC RADIATION (GLM)
!
      INTEGER IRAY_ATM
!
!     LAW OF TRACERS DEGRADATION
!
      INTEGER, ALLOCATABLE :: LOITRAC(:)
!
!     COEFFICIENT 1 FOR LAW OF TRACERS DEGRADATION
!     (1 IN CASE OF FUTURE LAW WITH MORE COEF.)
      DOUBLE PRECISION, ALLOCATABLE :: COEF1TRAC(:)
!
!     AED2 COUPLING
!
      INTEGER NWQVARS,NWQBEN,NWQDIAGS
!
!-----------------------------------------------------------------------
!
!     WATER QUALITY INDICES
!     (EVEN IF SOME ARE DUPLICATED IN TELEMAC)
!
!     TEMPERATURE
      INTEGER IND_T
!     SALINITY
      INTEGER IND_S
!     DISSOLVED O2
      INTEGER IND_O2
!     ORGANIC LOAD
      INTEGER IND_OL
!     NH4 LOAD
      INTEGER IND_NH4
!     PHYTO BIOMASS
      INTEGER IND_PHY
!     DISSOLVED PO4
      INTEGER IND_PO4
!     POR NON ASSIM
      INTEGER IND_POR
!     DISSOLVED NO3
      INTEGER IND_NO3
!     NOR NON ASSIM
      INTEGER IND_NOR
!     SUSPENDED LOAD
      INTEGER IND_SS
!     BED SEDIMENTS
      INTEGER IND_SF
!     MICRO POLLUTANT
      INTEGER IND_C
!     ABS. SUSP. LOAD.
      INTEGER IND_CSS
!     ABSORB. BED SED.
      INTEGER IND_CSF
!     THIS COULD ACTUALLY BE AN ARRAY
      INTEGER IND_AED2(MAXWQVAR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
! PART ADDED BY SERGIO
!     ZOO BIOMASS
! THIS IS THE IDENTIFIER OF ZOOPLANKTON INTO THE TRACERS BLOCK, TELEMAC
! ASSIGN IT A VALUE AUTOMATICALLY BY CALLING THE SUBROUTINE ADDTRACER
! INTO THE SUBROUTINE NAMETRAC_WAQTEL.F
      INTEGER IND_ZOO      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!
!-----------------------------------------------------------------------
!
!       3) LOGICAL
!
!-----------------------------------------------------------------------
!
!     IF YES, MASS-BALANCE
!
      LOGICAL WQBILMAS
!
!     IF YES, VALIDATION
!
      LOGICAL WQVALID
!
!     STEADY HYDRODYNAMICS
!
      LOGICAL  PERMA
!
!     IF TEMPERATURE IS VARIABLE
!
      LOGICAL YATEMP
!
!     AED2 COUPLING
!
      LOGICAL :: DEJA_SW = .FALSE.
      LOGICAL :: DEJA_CA = .FALSE.
!-----------------------------------------------------------------------
!
!       4) REALS
!
!-----------------------------------------------------------------------
!
!     WATER DENSITY
!
      DOUBLE PRECISION RO0 ! (RHO ZERO)
!
!     KINEMATIC VISCOSITY
!
      DOUBLE PRECISION VCE

!     LONGITUDINAL DISPERSION
!
      DOUBLE PRECISION LDISP
!
!     TRANSVERSAL DISPERSION
!
      DOUBLE PRECISION TDISP
!
!     WATER QUALITY VARIABLES PCO2
!
      DOUBLE PRECISION PCO2
!
!     WATER QUALITY VARIABLES PVAP
!
      DOUBLE PRECISION PVAP
!
!     WATER QUALITY VARIABLES RAY3
!
      DOUBLE PRECISION RAY3
!
!     SOLAR RADIATION FOR AED2
!
      TYPE(BIEF_OBJ), TARGET :: RAYAED2
!
!     WATER QUALITY VARIABLES NEBU
!
      DOUBLE PRECISION NEBU
!
!     WATER QUALITY VARIABLES - RELATIVE HUMIDITY
!
      DOUBLE PRECISION HREL
!
!     WATER QUALITY VARIABLES - RAINFALL
!
      DOUBLE PRECISION RAINFALL
!
!     WATER QUALITY VARIABLES - AIR TEMPERATURE
!
      TYPE(BIEF_OBJ), TARGET :: TAIR
!
!     WATER QUALITY VARIABLES - ???
!
      TYPE(BIEF_OBJ), TARGET :: TDEW
!
!     WATER QUALITY VARIABLES - ???
!
      TYPE(BIEF_OBJ), TARGET :: VISBI
!
!     VALUE OF THE AIR TEMPERATURE
!
      DOUBLE PRECISION TAIR_VALUE
!
!     WATER QUALITY VARIABLES ZSD (SECCHI LENGTH)
!
      DOUBLE PRECISION ZSD
!
!     WATER QUALITY VARIABLES NWIND (IF CONSTANT WIND IN SPACE)
!
      DOUBLE PRECISION NWIND
!
!     WATER QUALITY VARIABLES C14_ATM
!
      DOUBLE PRECISION C14_ATM
!
!     WATER QUALITY VARIABLES HTO_ATM
!
      DOUBLE PRECISION HTO_ATM
!
!     WEIR COEFFICIENT OF REAERATION
!
      DOUBLE PRECISION RSW
!
!     COEFFICIENT TO CALIBRATE THE ATMOSPHERE-WATER EXCHANGE MODEL
!
      DOUBLE PRECISION C_ATMOS
!
!     WATER QUALITY VARIABLE: EVAPORATION
!
      DOUBLE PRECISION EVAPORATION
!
! Water quality specific key-words
!
!
!     EUTRO PROCESS
!
!
!     CONSTANT OF DEGRADATION OF ORGANIC LOAD K120
!
      DOUBLE PRECISION K120
!
!     CONSTANT OF NITRIFICATION KINETIC K520
!
      DOUBLE PRECISION K520
!
!     O2 PRODUCED BY PHOTOSYNTHESIS (F FOR TRACER)
!
      DOUBLE PRECISION O2PHOTO
!
!     O2 CONSUMED BY NITRIFICATION (N FOR TRACER)
!
      DOUBLE PRECISION O2NITRI
!
!     BENTHIC DEMAND
!
      DOUBLE PRECISION DEMBEN
!
!      COEFFICIENT A AND B USED IN RS FORMULA
!
      DOUBLE PRECISION ABRS(2)
!
!     O2 SATURATION DENSITY OF WATER (CS)
!
      DOUBLE PRECISION O2SATU
!
!     SEDIMENTATION VELOCITY OF ORGANIC PHOSPHORUS
!
      DOUBLE PRECISION WPOR
!
!     SEDIMENTATION VELOCITY OF ORGANIC NITROGEN
!
      DOUBLE PRECISION WNOR
!
!     MAXIMUM ALGAL GROWTH
!
      DOUBLE PRECISION CMAX
!
!     COEF VEGETAL TURBIDITY WITHOUT PHYTOPLANKTON
!
      DOUBLE PRECISION KPE
!
!     PARAMETER OF CALIBRATION OF SMITH FORMULA
!
      DOUBLE PRECISION IK
!
!     CONSTANT OF HALF-SATURATION WITH PHOSPHATE
!
      DOUBLE PRECISION KP
!
!     CONSTANT OF HALF-SATURATION WITH PHOSPHA
!
      DOUBLE PRECISION KN
!
!     ALGAL COEFF OF TOXICITY (ALPHA; ALPHA2)
!
      DOUBLE PRECISION CTOXIC(2)
!
!     RESPIRATION RATE OF ALGAL BIOMASS (RP)
!
      DOUBLE PRECISION TRESPIR
!
!     PROPORTION OF PHOSPHORUS WITHIN PHYTO CELLS (FP)
!
      DOUBLE PRECISION PROPHOC
!
!     PERCENTAGE OF PHYSPHORUS ASSIMILABLE IN DEAD PHYTO
!
      DOUBLE PRECISION DTP
!
!     RATE OF TRANSFORMATION OF POR INTO PO
!
      DOUBLE PRECISION K320
!
!     PROPORTION OF NITROGEN WITHIN PHYTO CELLS (FN)
!
      DOUBLE PRECISION PRONITC
!
!     PERCENTAGE OF NITROGEN ASSIMILABLE IN DEAD PHYTO (DTN)
!
      DOUBLE PRECISION PERNITS
!
!     RATE OF TRANSFORMATION OF NOR INTO NO3
!
      DOUBLE PRECISION K360
!
!     COEF OF ALGAL MORTALITY (M1 AND M2)
!
      DOUBLE PRECISION CMORALG(2)
!
!     SEDIMENTATION VELOCITY OF ORGANIC LOAD
!
      DOUBLE PRECISION WLOR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
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
!     FEEDING RATE OF THE ZOOPLANKTON FROM PHYTOPLANKTON, ZOOPLANKTON AN
!     D ORGANIC LOAD
      TYPE(BIEF_OBJ), TARGET :: ZOOGS1
      TYPE(BIEF_OBJ), TARGET :: ZOOGS2
      TYPE(BIEF_OBJ), TARGET :: ZOOGS3
!
!     ASSIMILATION FEEDING OF ZOOPLANKTON FROM PHYTOPLANKTON, ZOOPLANKTO
!     N AND ORGANIC LOAD
      DOUBLE PRECISION, DIMENSION(3) :: ZOOOME
!
!     ZOOPLANKTON DIET FROM PHYTOPLANKTON, ZOOPLANKTON AND
!     ORGANIC LOAD
      DOUBLE PRECISION, DIMENSION(3) :: ZOORHOS
!
!     WORKING ARRAYS FOR DEAL WITH ZOOPLANKTON
      TYPE(BIEF_OBJ), TARGET :: ZOOT1
      TYPE(BIEF_OBJ), TARGET :: ZOOT2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!
!     BIOMASS PROCESS
!
!
!
!     EFFECTS OF PHOSPHORIOUS AND NITROGENIOUS NUTRIMENTS ON ALGAE GROWT
!
!      DOUBLE PRECISION LNUT
!
!     DENSITY OF SUNSHINE FLUX ON WATER SURFACE
!
      DOUBLE PRECISION I0
!
!
!     O2 PROCESS
!
!
!     CONSTANT OF DEGRADATION OF ORGANIC LOAD K1
!
      DOUBLE PRECISION K1
!
!     COEFFICIENT OF REAERATION K2 (if constant)
!
      DOUBLE PRECISION K22
!
!     CONSTANT OF NITRIFICATION KINETIC K4
!
      DOUBLE PRECISION K44
!
!     PHOTOSYNTHESIS P
!
      DOUBLE PRECISION PHOTO
!
!     VEGATAL RESPIRATION
!
      DOUBLE PRECISION  RESP
!
!     WATER TEMPERATURE
!
      DOUBLE PRECISION WATTEMP
!
!
!     MICROPOL PROCESS (MOSTLY COMMON WITH EUTRO)
!
!
!     EROSION RATE
!
      DOUBLE PRECISION ERO
!
!     SEDIMENTATION CRITICAL STRESS
!
      DOUBLE PRECISION TAUS
!
!     CRITICAL STRESS OF RESUSPENSION
!
      DOUBLE PRECISION TAUR
!
!     SEDIMENT SETTLING VELOCITY
!
      DOUBLE PRECISION VITCHU
!
!     COEFF OF DISTRIBUTION (KD)
!
      DOUBLE PRECISION CDISTRIB
!
!     CONSTANT OF DESORPTION KINETIC (KDESORP)
!
      DOUBLE PRECISION KDESORP
!
!
!
!     THERMIC PROCESS (COMMON VARIABLES WITH EWCHANGE_WITH_ATMOSPHERE MODULE)
!
!
!     WATER SPECIFIC HEAT
!
      DOUBLE PRECISION CP_EAU
!
!     AIR SPECIFIC HEAT
!
      DOUBLE PRECISION CP_AIR
!
!     COEFF OF CLOUDING RATE
!
      DOUBLE PRECISION COEF_K
!
!     COEFFICIENTS FOR CALIBRATING ATMOSPHERIC RADIATION
!
      DOUBLE PRECISION EMA
!
!     COEFFICIENTS FOR CALIBRATING FREE SURFACE RADIATION
!
      DOUBLE PRECISION EMI_EAU
!
!     BOLTZMANN CONSTANT (wM-2K-4)
!
      DOUBLE PRECISION, PARAMETER :: BOLTZ=5.67D-8
!
!     CONVERSION SECOND TO DAYS
!
      DOUBLE PRECISION, PARAMETER :: SECTODAY=1.D0/86400.D0
!
!     COEFFICIENTS OF RAERATION FORMULA
!
      DOUBLE PRECISION CFAER(2)
!
!
!     MES PROCESS
!
!
!     EXPONENTIAL DESINTEGRATION CONSTANT (LAMBD)
!
      DOUBLE PRECISION CCSEDIM
!
!     EVAPORATION
!
      DOUBLE PRECISION EVAPOR
!
!
!     AED2 PROCESS
!
!
      DOUBLE PRECISION, ALLOCATABLE :: EXTCAED2(:,:),FLUXAED2(:,:,:)
!
!-----------------------------------------------------------------------
!
!       5) STRINGS
!
!-----------------------------------------------------------------------
!
!     TITLE OF STUDY
!
      CHARACTER(LEN=72) TITWAQCAS
!
!     COPY OF SUBMIT STRINGS IN THE DICTIONARY
!
      CHARACTER(LEN=PATH_LEN) SUBMIT(4,300)

!     MAXIMUM OF LOGICAL UNITS NUMBERS
!
      INTEGER, PARAMETER :: MAXLU_WAQ = 14
!
!     BIEF_FILES STRUCTURES
!
      TYPE(BIEF_FILE) :: WAQ_FILES(MAXLU_WAQ)
!
      SAVE
!
      END MODULE DECLARATIONS_WAQTEL
                        
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
     &  IND_PHY,IND_PO4,IND_POR,IND_NO3,IND_NOR,IND_NH4,IND_OL,IND_O2,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDED PART BY SERGIO
! INCLUSING ZOOPLANKTON RELATED VARIABLES
     & IND_ZOO,MUZ,GAMMAZ,ZOOBK,ZOOG,BETAPC,BETAM3L,BETANC,BETAO2C,
     & ZOOGS1,ZOOGS2,ZOOGS3,ZOOOME,ZOORHOS,ZOOT1,ZOOT2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
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
      INTEGER                     :: I,II,JJ
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
!      WRITE(LU,*)'IND_PHY','IND_O2','IND_ZOO',IND_PHY,IND_O2,IND_ZOO
!      CALL PLANTE(1)
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
