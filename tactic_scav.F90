MODULE M_TACTIC_SCAV
 
    CONTAINS

    SUBROUTINE TACTIC_SCAV &
       ! -- input
        & (  YDCST_ACLIB, YPAER_ACLIB, YAEROSOLSK, &
        &   KIDIA   , KFDIA   , KLON    , KLEV    , &
        &   KLUOUT, & 
        &   PTSPHY  , &
        &   KWDEPIN, KWDEPOPT1, &
        &   PDP     , PTP     , &
        &   PDZ     , PRHO    , &
        &   PFLXCR  , PFLXCS  , PFLXSR  , PFLXSS  , &
       ! &   PCLWATCL, PCLWATCS, PCLWATSL, PCLWATSS, &
       ! &   PCLCOVC , PCLCOVS , PCLCOV  , &
        &   PCLCOV  , PCLWATL , PCLWATS , &
        &   KAER    , PAERO   , PTAERI  , &
       ! -- output
        &   PTAERO  , PFAERO , &
       ! &  ZPTAEROBCR, ZPTAEROBCS, ZPTAEROINR, ZPTAEROINS, ZPTAERORER , ZPTAERORES &
        &  PFAEROBCR, PFAEROBCS, PFAEROINR, PFAEROINS, PFAERORER , PFAERORES &
        & )
    
    
    
    !*** * TACTIC_SCAV* - IN-CLOUD SCAVENGING OF TRACERS
    
    
    
    !**   INTERFACE.
    !     ----------
    !          *TACTIC_SCAV* IS CALLED FROM *APLPAR*.
    
    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*
    !        FROM O. BOUCHER (LOA, 1998-03) in_scav
    
    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 2004-05-07
    !        P.BECHTOLD 2009-04-22 Clean/optimise and rewrite to externalise
    !                              scavenging coefficient
    !     End Modifications - CY37t1_bf.02
    !
    !        M.Michou 2012-02-10   Prognostic aerosols : modified for performance on a vectorial computer
    !        M.Michou 2014-03      Some CY38R2_MACC modifs reported here. Use on a scalar computer
    !        P.Nabat  2016-04      Optimization for the use of interactive aerosols
    !        J.Bock   2018-06      Major update, merge in-cloud and below-cloud scavenging routines
    !-----------------------------------------------------------------------
    
    USE PARKIND1  ,ONLY : JPIM  , JPRB
    USE YOMHOOK   ,ONLY : LHOOK , DR_HOOK
    USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
    USE MODD_PARAM_AER_ACLIB, ONLY :  TPARAM_AER_ACLIB, TAEROSOL
    
    IMPLICIT NONE
    
    !-----------------------------------------------------------------------
    
    !*       0.1  ARGUMENTS
    !             ---------
    TYPE(TCST_ACLIB),       INTENT(IN) :: YDCST_ACLIB
    TYPE(TPARAM_AER_ACLIB), INTENT(IN) :: YPAER_ACLIB
    TYPE(TAEROSOL),         INTENT(IN) :: YAEROSOLSK
    INTEGER(KIND=JPIM),     INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV
    INTEGER(KIND=JPIM),     INTENT(IN) :: KLUOUT 
    INTEGER(KIND=JPIM),     INTENT(IN) :: KAER
    INTEGER(KIND=JPIM),     INTENT(IN) :: KWDEPIN
    INTEGER(KIND=JPIM),     INTENT(IN) :: KWDEPOPT1
    
    REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY
    REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV)     ! Delta pressure [Pa]
    REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)     ! Temperature [K]
    REAL(KIND=JPRB),INTENT(IN)    :: PDZ(KLON,KLEV)     ! Delta altitude [m]
    REAL(KIND=JPRB),INTENT(IN)    :: PRHO(KLON,KLEV)    ! Dry air density [kg.m-3]
    
    REAL(KIND=JPRB),INTENT(IN)    :: PFLXCR(KLON,0:KLEV), PFLXCS(KLON,0:KLEV) ! Convective precip (rain/snow)
    REAL(KIND=JPRB),INTENT(IN)    :: PFLXSR(KLON,0:KLEV), PFLXSS(KLON,0:KLEV) ! Stratiform precip (rain/snow)
    REAL(KIND=JPRB),INTENT(IN)    :: PCLCOV(KLON,KLEV) , PCLWATL(KLON,KLEV) , PCLWATS(KLON,KLEV)
    REAL(KIND=JPRB),INTENT(IN)    :: PAERO(KLON,KLEV)  , PTAERI(KLON,KLEV)
    
    REAL(KIND=JPRB),INTENT(OUT) :: PFAERO(KLON), PTAERO(KLON,KLEV)
    REAL(KIND=JPRB),INTENT(OUT) :: PFAEROBCR(KLON), PFAEROBCS(KLON)  ! BC scav tendencies due to rain, snow vertical integral
    REAL(KIND=JPRB),INTENT(OUT) :: PFAEROINR(KLON), PFAEROINS(KLON)  ! IN scav tendencies due to rain, snow vertical integral
    REAL(KIND=JPRB),INTENT(OUT) :: PFAERORER(KLON), PFAERORES(KLON)  ! reevap scav tendencies due to rain, snow vertical integral
    
    !*       0.5   LOCAL VARIABLES
    !              ---------------
    
    INTEGER(KIND=JPIM) :: JK, JL
    
    REAL(KIND=JPRB) :: ZAERI(KLON,KLEV)   ! aerosol mass mixing ratio         [kg.kg-1]
    REAL(KIND=JPRB) :: ZPCLWAT(KLON,KLEV) ! total cloud water (liquid+solid)
    
    REAL(KIND=JPRB) :: ZPFLXR(KLON,0:KLEV) ! total rain flux (strati+conv)    [kg.m-2.s-1]
    REAL(KIND=JPRB) :: ZPFLXS(KLON,0:KLEV) ! total snow flux (strati+conv)    [kg.m-2.s-1]
    
    ! actual temperature of the parameterised type of cloud, depending on parameterisation choice
    !REAL(KIND=JPRB) :: ZTP_LCLDCV(KLON,KLEV), ZTP_MCLDCV(KLON,KLEV), ZTP_SCLDCV(KLON,KLEV)
    REAL(KIND=JPRB) :: ZTP_LCLD(KLON,KLEV), ZTP_MCLD(KLON,KLEV), ZTP_SCLD(KLON,KLEV)
    
    ! diagnose the presence (1) or absence (0) of liq/mix/sol cloud, depending on parameterisation choice
    !REAL(KIND=JPRB) :: ZIS_LCLDCV(KLON,KLEV), ZIS_MCLDCV(KLON,KLEV), ZIS_SCLDCV(KLON,KLEV)
    REAL(KIND=JPRB) :: ZIS_LCLD(KLON,KLEV), ZIS_MCLD(KLON,KLEV), ZIS_SCLD(KLON,KLEV)
    
    REAL(KIND=JPRB) :: ZALPHAR, ZALPHAS
    !REAL(KIND=JPRB) :: ZBETA, ZBETAI, ZBETAK, ZBETAR
    REAL(KIND=JPRB) :: ZBETARAIN, ZBETAIRAIN, ZBETAKRAIN, ZBETARRAIN
    REAL(KIND=JPRB) :: ZBETASNOW, ZBETAISNOW, ZBETAKSNOW, ZBETARSNOW
    REAL(KIND=JPRB) :: ZC1, ZC2
    REAL(KIND=JPRB) :: ZDP, ZDPORG ! Delta Pres, Delta Pres over RG
    REAL(KIND=JPRB) :: ZDX, ZDXR, ZDXS
    REAL(KIND=JPRB) :: ZEPS, ZF
    REAL(KIND=JPRB) :: ZFUNCR, ZFUNCS
    REAL(KIND=JPRB) :: ZINCRR, ZINCRS
    REAL(KIND=JPRB) :: ZSWATFRAC
    REAL(KIND=JPRB) :: ZEPSQLIQ, ZEPSFLX, Z1TSPHY, Z1RG, ZSCAV
    REAL(KIND=JPRB) :: ZPR, ZPS ! rain/snow flux at full level (average of upper and lower half levels) [kg(water).m-2.s-1]
    
    REAL(KIND=JPRB) :: ZPTAEROBCR(KLON,KLEV), ZPTAEROBCS(KLON,KLEV)             ! BC scav tendencies due to rain, snow
    REAL(KIND=JPRB) :: ZPTAEROINR(KLON,KLEV), ZPTAEROINS(KLON,KLEV)             ! IN scav tendencies due to rain, snow
    REAL(KIND=JPRB) :: ZPTAERORER(KLON,KLEV), ZPTAERORES(KLON,KLEV) ! reevap scav tendencies due to rain, snow
    REAL(KIND=JPRB) :: ZPFAEROR(KLON),        ZPFAEROS(KLON)        ! total scav tendencies due to rain, snow vertical integral
    
    ! scavenging coefficients for liquid/mixed/solid clouds
    REAL(KIND=JPRB) :: ZSCAV_L, ZSCAV_M, ZSCAV_S
    REAL(KIND=JPRB) :: ZMINCLCOV, ZCLCOV
    
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_SCAV:TACTIC_SCAV',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
    
    !- PAERO  in kg kg-1
    !  PFLXx  in kg m-2 s-1
    !  PCLWAT  in kg kg-1
    !  ZSCAV  a fraction N.D.
    !  PTAERI in kg m-3 s-1
    !  PTAERO in KG m-3 s-1
    
    ! NB: PCLWAT is the in-cloud water mixing ratio
    
    Z1RG=1.0_JPRB/YDCST_ACLIB%RG 
    Z1TSPHY=1.0_JPRB/PTSPHY
    
    ZEPSFLX =1.E-18_JPRB
    ZEPSQLIQ=1.E-18_JPRB
    
    !-- minimum cloud cover (same as in ACNEBN: ZEPSNEB)
    ZMINCLCOV=1.E-10_JPRB
    
    ! BC scav
    ZALPHAR = YAEROSOLSK%RSCAVBCR
    ZALPHAS = YAEROSOLSK%RSCAVBCS
    
    ZC1 = ZALPHAR/YPAER_ACLIB%RRAD_R_DROP/YDCST_ACLIB%RRHOLW
    ZC2 = ZALPHAS/YPAER_ACLIB%RRAD_S_CRYS/YDCST_ACLIB%RRHOIW
    
    
    ! initialise arrays
    PFAERO(KIDIA:KFDIA)     = 0._JPRB
    PTAERO(KIDIA:KFDIA,:)   = 0._JPRB
    ZPTAEROBCR(KIDIA:KFDIA,:) = 0._JPRB
    ZPTAEROBCS(KIDIA:KFDIA,:) = 0._JPRB
    ZPTAEROINR(KIDIA:KFDIA,:) = 0._JPRB
    ZPTAEROINS(KIDIA:KFDIA,:) = 0._JPRB
    ZPTAERORER(KIDIA:KFDIA,:) = 0._JPRB
    ZPTAERORES(KIDIA:KFDIA,:) = 0._JPRB
    PFAEROBCR(KIDIA:KFDIA) = 0._JPRB
    PFAEROBCS(KIDIA:KFDIA) = 0._JPRB
    PFAEROINR(KIDIA:KFDIA) = 0._JPRB
    PFAEROINS(KIDIA:KFDIA) = 0._JPRB
    PFAERORER(KIDIA:KFDIA) = 0._JPRB
    PFAERORES(KIDIA:KFDIA) = 0._JPRB
    ZPFAEROR(KIDIA:KFDIA)  = 0._JPRB
    ZPFAEROS(KIDIA:KFDIA)  = 0._JPRB
    
    ZSCAV = YAEROSOLSK%RSCAVIN
    ZSCAV_L = YAEROSOLSK%RSCAVINL
    ZSCAV_M = YAEROSOLSK%RSCAVINM
    ZSCAV_S = YAEROSOLSK%RSCAVINS
    
    
    
    ! Local arrays for precipitation (rain / snow), merge convective and stratiform arrays
    ZPFLXR(KIDIA:KFDIA,:) = MAX(0.0_JPRB,PFLXCR(KIDIA:KFDIA,:) + PFLXSR(KIDIA:KFDIA,:))
    ZPFLXS(KIDIA:KFDIA,:) = MAX(0.0_JPRB,PFLXCS(KIDIA:KFDIA,:) + PFLXSS(KIDIA:KFDIA,:))
    
    ! Local array for total cloud water (liquid + solid)
    ZPCLWAT(KIDIA:KFDIA,:) = MAX(0.0_JPRB,PCLWATL(KIDIA:KFDIA,:) + PCLWATS(KIDIA:KFDIA,:))
    
    ! Initial concentration
    ZAERI(KIDIA:KFDIA,1:KLEV) = PAERO(KIDIA:KFDIA,1:KLEV) + PTSPHY * PTAERI(KIDIA:KFDIA,1:KLEV)
    
    
    DO JK=1,KLEV
    
      !ZAERI(KIDIA:KFDIA,JK) = PAERO(KIDIA:KFDIA,JK) + PTSPHY * PTAERI(KIDIA:KFDIA,JK)
    
      DO JL=KIDIA,KFDIA
    
        ZDP =PDP(JL,JK)                           ! in Pa (kg m-1 s-2)
        ! jjb test a conversion using dry air
        !ZDPORG = ZDP*Z1RG
        ZDPORG = PRHO(JL,JK)*PDZ(JL,JK)
    
    ! ****************************************
    !  SECURITY CHECKS
    ! ****************************************
    
    
        ZBETAIRAIN=ZPFLXR(JL,JK)-ZPFLXR(JL,JK-1)   ! in kg m-2 s-1 (w/o cloud fraction effect)
        ZBETAISNOW=ZPFLXS(JL,JK)-ZPFLXS(JL,JK-1)   ! in kg m-2 s-1 (w/o cloud fraction effect)
    
    ! --------------------------------------------------------------------------------------------------------
    !  Re-evaporation
    ! --------------------------------------------------------------------------------------------------------
        IF(YPAER_ACLIB%LAER_REEVAP .AND. JK > 1)THEN
    ! RAIN part
    !    ZBETAI > 0 : scavenging.
    !    ZBETAI < 0 : re-evaporation. ZPFLX (top) must be non zero.
    
          IF (ZBETAIRAIN < 0._JPRB .AND. ZPFLXR(JL,JK-1) > ZEPSFLX) THEN
            ZBETARRAIN = ZBETAIRAIN / ZPFLXR(JL,JK-1)                   ! -- <0
          
            ZF = -ZBETARRAIN                                            !    >0
            ZEPS = (1. - EXP(-2.*ZF**0.5) * &
                   & (1. + 2.*ZF**0.5 + 2.*ZF + 4./3.*ZF**1.5)) &
                   & * (1. - ZF) + ZF**2
    
            IF(ZF < 0._JPRB .OR. ZF > 1._JPRB) THEN
              WRITE(KLUOUT,'("TACTIC_SCAV: problem with zf(1)",3I3,2E12.5)') JL,JK, KAER,ZF,ZEPS
              ZF=MIN(MAX(0._JPRB,ZF),1._JPRB)
            ENDIF
            IF(ZEPS < 0._JPRB .OR. ZEPS > 1._JPRB) THEN
              WRITE(KLUOUT,'("TACTIC_SCAV: problem with zeps(1)",3I3,2E12.5)') JL,JK, KAER,ZF,ZEPS
              ZEPS=MIN(MAX(0._JPRB,ZEPS),1._JPRB)
            ENDIF
    
            ! jjb test a conversion using dry air
            !ZDX = ZF*ZEPS * ZPFAEROR(JL) * RG/ZDP                     ! in kg kg-1 s-1  ! >0
            ZDX = ZF*ZEPS * ZPFAEROR(JL) / ZDPORG                     ! in kg kg-1 s-1  ! >0
            ZPTAERORER(JL,JK) = ZDX                                   ! in kg kg-1 s-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
            PFAERORER(JL) =  PFAERORER(JL) + ZDX*ZDPORG ! = PFAERORER(JL) + ZF*ZEPS*ZPFAEROR(JL)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
            ZPFAEROR(JL) =  (1._JPRB-ZF*ZEPS) * ZPFAEROR(JL)                            ! >0
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
            ! update aerosol concentration
            ZAERI(JL,JK) = ZAERI(JL,JK) + PTSPHY * ZPTAERORER(JL,JK)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
    
          ELSE ! ZBETAIRAIN >0
            ZPTAERORER(JL,JK) = 0._JPRB
          ENDIF
    
    ! SNOW part
    !    ZBETAI > 0 : scavenging.
    !    ZBETAI < 0 : re-evaporation. ZPFLX (top) must be non zero.
    
          IF (ZBETAISNOW < 0._JPRB .AND. ZPFLXS(JL,JK-1) > ZEPSFLX) THEN
            ZBETARSNOW = ZBETAISNOW / ZPFLXS(JL,JK-1)                   ! -- <0
    
            SELECT CASE(KWDEPOPT1)
            CASE (1)
              IF(ZBETARSNOW > 0._JPRB .OR. ZBETARSNOW < -1._JPRB) THEN
                WRITE(KLUOUT,'("TACTIC_SCAV: problem with zbetarsnow",3I3,2E12.5)') JL,JK, KAER,ZBETARSNOW,ZBETAISNOW
                ZBETARSNOW=MIN(MAX(-1._JPRB,ZBETARSNOW),0._JPRB)
              ENDIF
    
              ! jjb test a conversion using dry air
              !ZDX = - ZBETARSNOW * ZPFAEROS(JL) * RG/ZDP                ! in kg kg-1 s-1
              ZDX = - ZBETARSNOW * ZPFAEROS(JL) / ZDPORG                ! in kg kg-1 s-1
    
              ZPTAERORES(JL,JK) = ZDX                                   ! in kg kg-1 s-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
              PFAERORES(JL) =  PFAERORES(JL) + ZDX*ZDPORG
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
              
              ZPFAEROS(JL) =  (1._JPRB+ZBETARSNOW) * ZPFAEROS(JL)  ! mind the + sign : ZBETARSNOW <0
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
            CASE(2)
              ZF = -ZBETARSNOW                                            !    >0
              ZEPS = (1. - EXP(-2.*ZF**0.5) * &
                     & (1. + 2.*ZF**0.5 + 2.*ZF + 4./3.*ZF**1.5)) &
                     & * (1. - ZF) + ZF**2
    
              IF(ZF < 0._JPRB .OR. ZF > 1._JPRB) THEN
                WRITE(KLUOUT,'("TACTIC_SCAV: problem with zf(2)",3I3,4E12.5)') JL,JK, KAER,ZF,ZEPS,ZPFLXS(JL,JK),ZPFLXS(JL,JK-1)
                ZF=MIN(MAX(0._JPRB,ZF),1._JPRB)
              ENDIF
              IF(ZEPS < 0._JPRB .OR. ZEPS > 1._JPRB) THEN
                WRITE(KLUOUT,'("TACTIC_SCAV: problem with zeps(2)",3I3,4E12.5)') JL,JK, KAER,ZF,ZEPS,ZPFLXS(JL,JK),ZPFLXS(JL,JK-1)
                ZEPS=MIN(MAX(0._JPRB,ZEPS),1._JPRB)
              ENDIF
    
              ! jjb test a conversion using dry air
              !ZDX = ZF*ZEPS * ZPFAEROS(JL) * RG/ZDP                     ! in kg kg-1 s-1  ! >0
              ZDX = ZF*ZEPS * ZPFAEROS(JL) / ZDPORG                     ! in kg kg-1 s-1  ! >0
    
              ZPTAERORES(JL,JK) = ZDX                                   ! in kg kg-1 s-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
              PFAERORES(JL) =  PFAERORES(JL) + ZDX*ZDPORG
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
              ZPFAEROS(JL) =  (1._JPRB-ZF*ZEPS) * ZPFAEROS(JL)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
            ENDSELECT
    
            ! update aerosol concentration
            ZAERI(JL,JK) = ZAERI(JL,JK) + PTSPHY * ZPTAERORES(JL,JK)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          ELSE ! ZBETAISNOW >0
            ZPTAERORES(JL,JK) = 0._JPRB
          ENDIF
    
        ELSE ! LAERREEVAP=.F. OR JK=1
          ZPTAERORER(JL,JK) = 0._JPRB
          ZPTAERORES(JL,JK) = 0._JPRB
    
        ENDIF ! LAERREEVAP
    
    ! --------------------------------------------------------------------------------------------------------
    !  Below cloud scavenging: only if precipitations are leaving the grid cell
    ! --------------------------------------------------------------------------------------------------------
    ! BC scavenging due to rain, update concentration afterwards
    ! -----------------------------------------------------------
        IF(ZPFLXR(JL,JK) > ZEPSFLX) THEN
          ! average rain flux at full level
          ZPR = 0.5 * (ZPFLXR(JL,JK-1)+ZPFLXR(JL,JK))
          ! leaving flux
          !ZPR = ZPFLXR(JL,JK)
    
          ZINCRR = 0.75 * ZPR*ZC1
    
          ZPTAEROBCR(JL,JK) = -ZINCRR * ZAERI(JL,JK)                 ! <0
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          PFAEROBCR(JL) = PFAEROBCR(JL) - ZPTAEROBCR(JL,JK)*ZDPORG   ! >0
          ZPFAEROR(JL)  = ZPFAEROR(JL)  - ZPTAEROBCR(JL,JK)*ZDPORG   ! >0
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
    ! update concentration
          ZAERI(JL,JK) = ZAERI(JL,JK) + PTSPHY*ZPTAEROBCR(JL,JK)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
        ELSE
          ZPTAEROBCR(JL,JK) = 0._JPRB
        ENDIF
    
    ! BC scavenging due to snow, update concentrations afterwards
    ! -----------------------------------------------------------
        IF(ZPFLXS(JL,JK) > ZEPSFLX) THEN
          ! average snow flux at full level
          ZPS = 0.5 * (ZPFLXS(JL,JK-1)+ZPFLXS(JL,JK))
          ! leaving flux
          !ZPS = ZPFLXS(JL,JK)
    
          ZINCRS = 0.75 * ZPS*ZC2
          IF (ISNAN(ZINCRS)) print*, "tactic_scav nan", JL, JK, ZINCRS, ZPS, ZC2
          ZPTAEROBCS(JL,JK) = -ZINCRS * ZAERI(JL,JK)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          PFAEROBCS(JL) = PFAEROBCS(JL) - ZPTAEROBCS(JL,JK)*ZDPORG
          ZPFAEROS(JL)  = ZPFAEROS(JL)  - ZPTAEROBCS(JL,JK)*ZDPORG
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
    ! update concentration
          ZAERI(JL,JK) = ZAERI(JL,JK) + PTSPHY*ZPTAEROBCS(JL,JK)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
        ELSE
          ZPTAEROBCS(JL,JK) = 0._JPRB
        ENDIF
    
    
    ! --------------------------------------------------------------------------------------------------------
    !  In cloud scavenging
    ! --------------------------------------------------------------------------------------------------------
    
        !-- This security is maybe useless...
        ZCLCOV = MIN(MAX( 0._JPRB, PCLCOV(JL,JK)), 1._JPRB)
    
        IF (ZPCLWAT(JL,JK) > 0._JPRB) THEN
          ZSWATFRAC = PCLWATS(JL,JK) / ZPCLWAT(JL,JK)
        ELSE
          ZSWATFRAC = 0._JPRB
        ENDIF
    
    ! RAIN part
    !    ZBETAI > 0 : scavenging. CLCOV and CLWAT must be positive (larger than threshold) to compute
    !    ZBETAI < 0 : re-evaporation (see next section)
    
        IF (ZBETAIRAIN > 0._JPRB .AND. ZPCLWAT(JL,JK) > ZEPSQLIQ) THEN
          ZBETAKRAIN = ZBETAIRAIN                                         ! in kg m-2 s-1
          ZBETARAIN  = ZBETAKRAIN / (ZPCLWAT(JL,JK)*ZDPORG)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
    !-- Apply the scavenging coefficient
          SELECT CASE(KWDEPIN)
          ! original version
          CASE (1)
            ZFUNCR=EXP(-ZSCAV*ZBETARAIN*PTSPHY)                           ! N.D.       (N.D.)
          ! option 2: use only the liquid phase coef
          CASE (2)
            ZFUNCR=EXP(-ZSCAV_L*ZBETARAIN*PTSPHY)
          ! option 3: use ice/mixed-phase/liquid coefficients depending on the LARGE-SCALE ONLY cldwat ratio
          CASE (3)
            IF(ZSWATFRAC < 0.1_JPRB) THEN       ! Less than 10% solid water => assume liquid
              ZFUNCR=EXP(-ZSCAV_L*ZBETARAIN*PTSPHY)
    
              ! Diagnose the presence and temperature of parameterised type of cloud
              ZTP_LCLD(JL,JK) = PTP(JL,JK)
              ZIS_LCLD(JL,JK) = 1._JPRB
              ! no mixed/solid clouds here
              ZTP_MCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZTP_SCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZIS_MCLD(JL,JK) = 0._JPRB
              ZIS_SCLD(JL,JK) = 0._JPRB
    
            ELSEIF(ZSWATFRAC < 0.9_JPRB) THEN   ! 10 - 90% solid water      => assume mixed
              ZFUNCR=EXP(-ZSCAV_M*ZBETARAIN*PTSPHY)
    
              ! Diagnose the presence and temperature of parameterised type of cloud
              ZTP_MCLD(JL,JK) = PTP(JL,JK)
              ZIS_MCLD(JL,JK) = 1._JPRB
              ! no liquid/solid clouds here
              ZTP_LCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZTP_SCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZIS_LCLD(JL,JK) = 0._JPRB
              ZIS_SCLD(JL,JK) = 0._JPRB
    
            ELSE                                    ! More than 90% solid water => assume ice
              ZFUNCR=EXP(-ZSCAV_S*ZBETARAIN*PTSPHY)
    
              ! Diagnose the presence and temperature of parameterised type of cloud
              ZTP_SCLD(JL,JK) = PTP(JL,JK)
              ZIS_SCLD(JL,JK) = 1._JPRB
              ! no liquid/mixed clouds here
              ZTP_LCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZTP_MCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZIS_LCLD(JL,JK) = 0._JPRB
              ZIS_MCLD(JL,JK) = 0._JPRB
    
            ENDIF
          ENDSELECT
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
    
          ZDXR = ZAERI(JL,JK)*(1.0 - ZFUNCR) * ZCLCOV          ! in kg kg-1       ! >0
          ZPTAEROINR(JL,JK) = -ZDXR*Z1TSPHY                        ! in kg kg-1 s-1   ! <0
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          PFAEROINR(JL) = PFAEROINR(JL) - ZPTAEROINR(JL,JK)*ZDPORG ! in kg m-2 s-1
          ZPFAEROR(JL)  = ZPFAEROR(JL)  - ZPTAEROINR(JL,JK)*ZDPORG ! in kg m-2 s-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          ! update aerosol specific concentration
          ZAERI(JL,JK) = ZAERI(JL,JK) + PTSPHY*ZPTAEROINR(JL,JK)   ! in kg kg-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
        ELSE
          ZPTAEROINR(JL,JK) = 0._JPRB
          ! diagnostic
          ZTP_LCLD(JL,JK) = YDCST_ACLIB%RUNDEF
          ZIS_LCLD(JL,JK) = 0._JPRB
          ZTP_MCLD(JL,JK) = YDCST_ACLIB%RUNDEF
          ZIS_MCLD(JL,JK) = 0._JPRB
          ZTP_SCLD(JL,JK) = YDCST_ACLIB%RUNDEF
          ZIS_SCLD(JL,JK) = 0._JPRB
        ENDIF   
    
    
    ! SNOW part
    !    ZBETAI > 0 : scavenging. CLCOV and CLWAT must be positive (larger than threshold) to compute
    !    ZBETAI < 0 : re-evaporation
    
        IF (ZBETAISNOW > 0._JPRB .AND. ZPCLWAT(JL,JK) > ZEPSQLIQ) THEN
          ZBETAKSNOW = ZBETAISNOW                                         ! in kg m-2 s-1
          ZBETASNOW  = ZBETAKSNOW / (ZPCLWAT(JL,JK)*ZDPORG)
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
    !-- Apply the scavenging coefficient
          SELECT CASE(KWDEPIN)
          ! original version
          CASE (1)
            ZFUNCS=EXP(-ZSCAV*ZBETASNOW*PTSPHY)                           ! N.D.       (N.D.)
          ! option 2: use only the liquid phase coef
          CASE (2)
            ZFUNCS=EXP(-ZSCAV_L*ZBETASNOW*PTSPHY)
          ! option 3: use ice/mixed-phase/liquid coefficients depending on the LARGE-SCALE ONLY cldwat ratio
          CASE (3)
            IF(ZSWATFRAC < 0.1_JPRB) THEN       ! Less than 10% solid water => assume liquid
              ZFUNCS=EXP(-ZSCAV_L*ZBETASNOW*PTSPHY)
    
              ! Diagnose the presence and temperature of parameterised type of cloud
              ZTP_LCLD(JL,JK) = PTP(JL,JK)
              ZIS_LCLD(JL,JK) = 1._JPRB
              ! no mixed/solid clouds here
              ZTP_MCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZTP_SCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZIS_MCLD(JL,JK) = 0._JPRB
              ZIS_SCLD(JL,JK) = 0._JPRB
    
            ELSEIF(ZSWATFRAC < 0.9_JPRB) THEN   ! 10 - 90% solid water      => assume mixed
              ZFUNCS=EXP(-ZSCAV_M*ZBETASNOW*PTSPHY)
    
              ! Diagnose the presence and temperature of parameterised type of cloud
              ZTP_MCLD(JL,JK) = PTP(JL,JK)
              ZIS_MCLD(JL,JK) = 1._JPRB
              ! no liquid/solid clouds here
              ZTP_LCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZTP_SCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZIS_LCLD(JL,JK) = 0._JPRB
              ZIS_SCLD(JL,JK) = 0._JPRB
    
            ELSE                                    ! More than 90% solid water => assume ice
              ZFUNCS=EXP(-ZSCAV_S*ZBETASNOW*PTSPHY)
    
              ! Diagnose the presence and temperature of parameterised type of cloud
              ZTP_SCLD(JL,JK) = PTP(JL,JK)
              ZIS_SCLD(JL,JK) = 1._JPRB
              ! no liquid/mixed clouds here
              ZTP_LCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZTP_MCLD(JL,JK) = YDCST_ACLIB%RUNDEF
              ZIS_LCLD(JL,JK) = 0._JPRB
              ZIS_MCLD(JL,JK) = 0._JPRB
    
            ENDIF
          ENDSELECT
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          ZDXS = ZAERI(JL,JK)*(1._JPRB - ZFUNCS) * ZCLCOV          ! in kg kg-1
          ZPTAEROINS(JL,JK) = -ZDXS*Z1TSPHY                        ! in kg kg-1 s-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          PFAEROINS(JL) = PFAEROINS(JL) - ZPTAEROINS(JL,JK)*ZDPORG ! in kg m-2 s-1
          ZPFAEROS(JL)  = ZPFAEROS(JL)  - ZPTAEROINS(JL,JK)*ZDPORG ! in kg m-2 s-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
          ! update aerosol specific concentration
          ZAERI(JL,JK) = ZAERI(JL,JK) + PTSPHY*ZPTAEROINS(JL,JK)   ! in kg kg-1
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
    
        ELSE
          ZPTAEROINS(JL,JK) = 0._JPRB
          ! diagnostic
          ZTP_LCLD(JL,JK) = YDCST_ACLIB%RUNDEF
          ZIS_LCLD(JL,JK) = 0._JPRB
          ZTP_MCLD(JL,JK) = YDCST_ACLIB%RUNDEF
          ZIS_MCLD(JL,JK) = 0._JPRB
          ZTP_SCLD(JL,JK) = YDCST_ACLIB%RUNDEF
          ZIS_SCLD(JL,JK) = 0._JPRB
        ENDIF   
    
    
    ! --------------------------------------------------------------------------------------------------------
    !  Total tendency and flux
    ! --------------------------------------------------------------------------------------------------------
        
        IF (YPAER_ACLIB%LAER_REEVAP) THEN
          PTAERO(JL,JK) = PTAERI(JL,JK) + &
                        & ZPTAEROBCR(JL,JK) + ZPTAEROBCS(JL,JK) + & ! <0
                        & ZPTAEROINR(JL,JK) + ZPTAEROINS(JL,JK) + & ! <0
                        & ZPTAERORER(JL,JK) + ZPTAERORES(JL,JK)     ! >0
        ELSE
          PTAERO(JL,JK) = PTAERI(JL,JK) + &
                        & ZPTAEROBCR(JL,JK) + ZPTAEROBCS(JL,JK) + & ! <0
                        & ZPTAEROINR(JL,JK) + ZPTAEROINS(JL,JK)     ! <0
        ENDIF
    
    ! ****************************************
    !  SECURITY CHECKS (debug)
    ! ****************************************
      ENDDO
    ENDDO
    
    PFAERO(KIDIA:KFDIA) = ZPFAEROR(KIDIA:KFDIA) + ZPFAEROS(KIDIA:KFDIA)
    
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_SCAV:TACTIC_SCAV',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
    END SUBROUTINE TACTIC_SCAV

END MODULE M_TACTIC_SCAV
