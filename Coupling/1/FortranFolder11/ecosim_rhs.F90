!                    **********************
                     SUBROUTINE ECOSIM_RHS &
!                    **********************
!
     (ECOSM2,EQPERB,FIRSTDERV,TN,HPROP,TEXP,IND_PHY,IND_ZOO, &
       IND_OL,IND_ECO,NPOIN)
!
!***********************************************************************
! TELEMAC2D   V8P1R0
!***********************************************************************
!
!brief    COMPUTES ECOSIM RIGHT HAND SIDE TERMS
!
!history  CASTIBLANCO-BALLESTEROS S. A.
!+        23/12/2020
!+        V8P1
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ECOSM2         |<--| M2 VARIABLE TO COMPUTE RISK_RATIOS
!| EQPERB         |<--| QperB VARIABLE TO COMPUTE RISK_RATIOS
!| FIRSTDERV      |<--| HERE IT COMES EQUALS .TRUE.
!| TN             |<--| TRACERS BLOCK
!|                |   | PHY, ZOO = microg/L   ; OL = mgO2/L
!| HPROP          |<--| WATER DEPTH
!| TEXP           |<--| RHS EXPLICIT TERMS
!| IND_PHY        |<--| POSITIONAL INDEX OF PHYTOPLANKTON IN TN AND TEXP
!| IND_ZOO        |<--| POSITIONAL INDEX OF ZOOPLANKTON IN TN AND TEXP
!| IND_OL         |<--| POSITIONAL INDEX OF ORGANIC LOAD IN TN AND TEXP
!| IND_ECO        |<--| POSITIONAL INDEXES OF ECOGRUOPS IN TN AND TEXP
!| NPOIN          |<--| NUMBER OF POINTS WITHIN MESH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
      USE DECLARATIONS_TELEMAC2D, ONLY: AT
      USE DECLARATIONS_WAQTEL, ONLY: DENAVG
      USE statevartypesecosim
      USE statevartypesecopath
!
      IMPLICIT NONE
!
!  INTENT VARIABLES
      LOGICAL, INTENT(IN) :: FIRSTDERV
      TYPE(BIEF_OBJ), INTENT(IN)    :: TN, HPROP
      TYPE(BIEF_OBJ), INTENT(INOUT) :: ECOSM2,EQPERB,TEXP
      INTEGER, INTENT(IN) :: IND_PHY,IND_ZOO,IND_OL,NPOIN
      INTEGER, INTENT(IN), DIMENSION(nvars) :: IND_ECO

!
!  IN-SUBROUTINE VARIABLES
      INTEGER :: I
      INTEGER :: K
      DOUBLE PRECISION, DIMENSION(nvars) :: dumQperB, dumM2, BIOM
      DOUBLE PRECISION :: MICL_KGM3 = 0.000001D0
      DOUBLE PRECISION :: MGL_KGM3 = 0.001D0
      DOUBLE PRECISION :: EPSIL = 1E-14
      DOUBLE PRECISION, DIMENSION(nvars) :: lossSt
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    DEALING WITH TIME IN ECOSIM
!
      IF(AT.GT.(SecondsPerMonth*imonth)) THEN
!
      !WRITE(LU,*) "UP",step, imonth, n, AT
!
      !step = step + 1      !ADVANCING STEP
      !time = time + tstep  !ADVANCING ECO-TIME
      !n = n + 1            !ADVANCING N INDEX
!
      ! Clean monthly stanza variables
      BBAvg(:)   = 0
      LossAvg(:) = 0
      EatenByAvg(:) = 0
      EatenOfAvg(:) = 0
      PredAvg(:)    = 0
!
      imonth = imonth + 1
      !n = 1
      !
      !Update Stanzas change made in time_ecosim (time section of ecosim)
      !
        UpdateStanzas = .TRUE.
        DENAVG = DENAVG + 1
        !WRITE(*,*) 'DENAAAAAAAAAAAAVG ', 'UP', DENAVG
      ELSE
      !
      !Update Stanzas change made in time_ecosim (time section of ecosim)
      !
        UpdateStanzas = .FALSE.
        DENAVG = DENAVG + 1
      END IF
!
      !IF(AT.GT.(SecondsPerMonth*(imonth-1)
     !&   + n*(SecondsPerMonth/StepsPerMonth))) THEN

      !WRITE(LU,*) "DOWN",step, imonth, n, AT

      !step = step + 1      !ADVANCING STEP
      !time = time + tstep  !ADVANCING ECO-TIME
      !n = n + 1            !ADVANCING N INDEX

      !ENDIF
      time = AT/(SecondsPerMonth*12)
!      WRITE(*,*) 'TIIIME ', time
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      ! CLEANING ECOSM2 AND EQPERB
      CALL OS('X=0     ',X=ECOSM2)
      CALL OS('X=0     ',X=EQPERB)
!
      ! GOING FOR DERIVS TO COMPUTE RHS TERMS
      DO K = 1,NPOIN

       dumQperB = 0D0
       dumM2 = 0D0
       BIOM = 0D0
       IF(HPROP%R(K).GT.EPSIL) THEN
       !
       ! STORING BIOMASS VALUES IN BIOM VECTOR
       !
        DO I = 1,nvars
         IF(I.EQ.nvars) THEN
           BIOM(I) = TN%ADR(IND_OL)%P%R(K)*HPROP%R(K)*MGL_KGM3
         ELSE IF(I.EQ.(nvars - 1)) THEN
           BIOM(I) = TN%ADR(IND_PHY)%P%R(K)*HPROP%R(K)*MICL_KGM3
         ELSE IF(I.EQ.vcols) THEN
           BIOM(I) = TN%ADR(IND_ZOO)%P%R(K)*HPROP%R(K)*MICL_KGM3
         ELSE
           BIOM(I) = TN%ADR(IND_ECO(I))%P%R(K)
         END IF

        END DO
       !
       ! CALL TO DERIVS
       !
        CALL DERIVS &
        (nvars, ndetritus, vrows, vcols, imonth, time, BIOM, xdot, &
        biomeq, loss, integrate, ep_data, ep_detfate, flow2detritus, &
        det_export_rate, es_data, es_vul, arena, NutFree, NutTot, &
        NutBiom, NutFreeBase, NutMin, detritus_no, boolFN, boolFPP, &
        FirstTime,dumQperB,dumM2,FIRSTDERV)
!
       ! WRITE(*,*) 'XDOOOOOT ', xdot
       !
       ! GETTING NEW TERMS TO TEXP
       !
        DO I = 1,nvars
          IF(I.EQ.nvars) THEN
            xdot(I) = xdot(I)/MGL_KGM3
          ELSE IF((I.EQ.vcols).OR.(I.EQ.(nvars-1))) THEN
            xdot(I) = xdot(I)/MICL_KGM3
          END IF
        END DO
!        WRITE(LU,*) 'XDOOOOOOT ', xdot
        DO I = 1,nvars
         IF(I.EQ.nvars) THEN
          TEXP%ADR(IND_OL)%P%R(K) = TEXP%ADR(IND_OL)%P%R(K) - xdot(I)
         ELSE IF(I.EQ.(nvars - 1)) THEN
          TEXP%ADR(IND_PHY)%P%R(K) = TEXP%ADR(IND_PHY)%P%R(K) - xdot(I)
         ELSE IF(I.EQ.vcols) THEN
          TEXP%ADR(IND_ZOO)%P%R(K) = TEXP%ADR(IND_ZOO)%P%R(K) - xdot(I)
         ELSE
         TEXP%ADR(IND_ECO(I))%P%R(K) = xdot(I)
         END IF
!    
        END DO
       END IF
      
!
!     ADDED FOR COMPATIBILITY WITH ECOSPACE-F
!
      do i = 1, nvars
        EatenByAvg(i)  = EatenByAvg(i) + es_data(i)%qq
        EatenOfAvg(i)  = EatenOfAvg(i) + es_data(i)%M2
        PredAvg(i)     = PredAvg(i) + es_data(i)%pred
        lossSt(i)      = loss(i)
      end do
! averaging here for multistanza calculations
! and updating foraging times
      do i = 1, nvars
    
        BBAvg(i)   = BBAvg(i) + BIOM(i)
        LossAvg(i) = LossAvg(i) + loss(i)
        
        if (UpdateStanzas .eqv. .true.) then
         BBAvg(i)        = BBAvg(i) / DENAVG
         LossAvg(i)      = LossAvg(i) / DENAVG
         EatenByAvg(i)   = EatenByAvg(i) / DENAVG
         EatenOfAvg(i)   = EatenOfAvg(i) / DENAVG
         PredAvg(i)      = PredAvg(i) / DENAVG
         lossSt(i)       = lossSt(i) / DENAVG
        end if
            !WRITE(*,*) 'BB ',B(i)
      end do
      
      if (UpdateStanzas .eqv. .true.) then
    
       if (nstanzas > 0) then
        ! below was lossSt but corrected upon EwE6 rk4 routine check
        !WRITE(*,*) 'MULTISTA ', 'UP', DENAVG
        call multistanza (nvars,nstanzas,BIOM,ms_data,es_data, &
          es_ms_data,BBAvg,LossAvg)
        !call multistanza (nvars,nstanzas,BBAvg,ms_data,es_data, &
        !  es_ms_data,BBAvg,LossAvg)
        !WRITE(*,*) 'MULTISTA ', 'DOWN', DENAVG
       end if
        
       ! Update foraging times at the end of each time step
       call updateForagingTimes (nvars,integrate,ep_data,es_data, &
         EatenByAvg,EatenOfAvg,BBAvg,PredAvg)
      
       ! SETTING TO ZERO DENAVG AGAIN
       !WRITE(*,*) 'DENAAAAAAAAAAAAVG ', 'DOWN1', DENAVG
       DENAVG = 0
       !WRITE(*,*) 'DENAAAAAAAAAAAAVG ', 'DOWN2', DENAVG
      end if
      
      xdot = 0D0
      loss = 0D0
      biomeq = 0D0
      lossSt = 0D0

      END DO

      RETURN
      END SUBROUTINE
