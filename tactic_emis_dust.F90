MODULE M_TACTIC_EMIS_DUST

    CONTAINS

    SUBROUTINE TACTIC_EMIS_DUST &
        !---inputs
        & (YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
        &  KIDIA, KFDIA, KLON, KAERO, KLUOUT, &
        &  KDDUST, LDUST_SFXSOIL, &
        &  PBARESOIL, PLSM , PSNS, PAERDEP, PAERSCC, PZ0M, PWIND, PRHOS, PWS1, &
        !---outputs
        &  PFLXDUST)

        !*** * AER_DUST_TACTIC* - SOURCE TERMS FOR DUST AEROSOLS
        
        !**   INTERFACE.
        !     ----------
        !          *AER_DUST_TACTIC* IS CALLED FROM *AER_SRC*.
        
        !     AUTHOR.
        !     -------
        !        Original version
        !        P. Nabat August 2023
        
        !     SOURCE.
        !     -------
        
        !     MODIFICATIONS.
        !-----------------------------------------------------------------------
        
        USE PARKIND1  , ONLY : JPIM, JPRB
        USE YOMHOOK   , ONLY : LHOOK, DR_HOOK

        USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
        USE MODD_PARAM_AER_ACLIB, ONLY : TAEROSOL, TPARAM_AER_ACLIB

        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        
        !*       0.1   ARGUMENTS
        !              ---------
        
        !---input fields
        ! STRUCTURES
        TYPE(TCST_ACLIB),                INTENT(IN)  :: YDCST_ACLIB
        TYPE(TPARAM_AER_ACLIB),          INTENT(IN)  :: YPAER_ACLIB
        TYPE(TAEROSOL),DIMENSION(KAERO), INTENT(IN)  :: YAEROSOLS

        INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
        INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
        INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLUOUT

        INTEGER(KIND=JPIM),INTENT(IN)    :: KDDUST
        LOGICAL,INTENT(IN)               :: LDUST_SFXSOIL
        
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PBARESOIL(KLON) ! =PFRTI(KLON,8)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNS(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERDEP(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERSCC(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHOS(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS1(KLON)
        
        REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLXDUST(KLON,KAERO)
        
        !-----------------------------------------------------------------------

        !-------------------------  Dust constants--------------
        INTEGER(KIND=JPIM), PARAMETER :: ISIZE_DD = 12
        INTEGER(KIND=JPIM), PARAMETER :: INATS_DD = 12
        INTEGER(KIND=JPIM), PARAMETER :: INSOIL_DD = 152
        INTEGER(KIND=JPIM), PARAMETER :: INSOILMODE_DD = 3
        REAL(KIND=JPRB),    PARAMETER :: ZZ0SAL = 3.E-05_JPRB   ! en m
        REAL(KIND=JPRB),    PARAMETER :: ZDS_KOK = 3.4_JPRB
        REAL(KIND=JPRB),    PARAMETER :: ZSIGMAS_KOK = 3.0_JPRB
        REAL(KIND=JPRB),    PARAMETER :: ZLAMBDA_KOK = 12.0_JPRB
        REAL(KIND=JPRB),    PARAMETER :: ZCV_KOK = 12.62_JPRB
        REAL(KIND=JPRB),    PARAMETER :: ZEPSDD = 1.0E-7_JPRB
        REAL(KIND=JPRB),    PARAMETER :: ZRVONKAR = 0.4_JPRB
        
        !   LOCAL VARIABLES
        !              ---------------
        INTEGER(KIND=JPIM) :: JAER, JL, JS , JN, JNT, JNM
        REAL(KIND=JPRB) :: ZSTOTAL , ZXK , ZXL , ZXM , ZXN, ZRWI, ZTOTV
        REAL(KIND=JPRB) :: ZALPHAPROP, ZUTMIN
        REAL(KIND=JPRB) :: ZRE, ZKK, ZDELDP, ZUTH, ZFDP1, ZFDP2, ZDDP
        REAL(KIND=JPRB) :: ZRWPWP, ZRWSAT
        REAL(KIND=JPRB) :: ZCLAY, ZSAND, ZSILT
        REAL(KIND=JPRB) :: ZSWETN, ZSWETN2
        REAL(KIND=JPRB) :: ZRHO_DD
        
        REAL(KIND=JPRB), DIMENSION(KLON) :: ZFW, ZFEFF, ZUSTAR, ZSWI1
        
        REAL(KIND=JPRB), DIMENSION(ISIZE_DD+1) :: ZAEROSIZE_DD

        REAL(KIND=JPRB), DIMENSION(INSOIL_DD) :: ZZ0S
        REAL(KIND=JPRB), DIMENSION(INATS_DD,INSOIL_DD) :: ZSREL
        REAL(KIND=JPRB), DIMENSION(INSOIL_DD) :: ZSS, ZDP_ARRAY
        REAL(KIND=JPRB), DIMENSION(INATS_DD,ISIZE_DD,KLON) :: ZRSFROWSUB
        REAL(KIND=JPRB), DIMENSION(INATS_DD,KLON,KAERO) :: ZRSFROWT
        REAL(KIND=JPRB), DIMENSION(KLON,INSOIL_DD) :: ZUSTARTS, ZUSTART,ZFEFFS
        REAL(KIND=JPRB), DIMENSION(ISIZE_DD) :: ZFRAC
        REAL(KIND=JPRB), DIMENSION(INATS_DD,KLON) :: ZFSOIL, ZFTEX
        REAL(KIND=JPRB), DIMENSION(KAERO+1):: ZTRSIZE_DD
        REAL(KIND=JPRB), DIMENSION(INSOILMODE_DD,INATS_DD) :: ZMMD_SOIL
        REAL(KIND=JPRB), DIMENSION(INSOILMODE_DD,INATS_DD) :: ZSIGMA_SOIL
        REAL(KIND=JPRB), DIMENSION(INSOILMODE_DD,INATS_DD) :: ZPCENT_SOIL
        
        LOGICAL :: LLDUST(KLON), LLDDACT(KLON)
        
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_EMIS_DUST:TACTIC_EMIS_DUST',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
        
        !*       0.1   INITIALIZATIONS
        !              --------------------------
        
        ZRWPWP=0.1715 ! valeur recalculée
        ZRWSAT=0.47225 ! valeur recalculée (voir main/sur/module/cptave_mod.F90)
        ZSWI1(KIDIA:KFDIA) = 0.  ! valeur bidon à remplacer par ce qui vient de SURFEX
        
        ! Initialization
        DO JAER=1,KAERO
           ZTRSIZE_DD(JAER) = YAEROSOLS(JAER)%RDIAMIN*1.E6  ! conversion in um
        ENDDO
        ZTRSIZE_DD(KAERO+1) = YAEROSOLS(KAERO)%RDIAMAX*1.E6  ! conversion in um


        ! Parameters of Marticorena scheme
        IF (KDDUST==2) THEN
          ZMMD_SOIL(:,1)  = (/ 690.0 , 210.0 , 10.0 /)
          ZMMD_SOIL(:,2)  = (/ 690.0 , 210.0 ,  0.0 /)
          ZMMD_SOIL(:,3)  = (/ 690.0 , 210.0 ,  0.0 /)
          ZMMD_SOIL(:,4)  = (/ 520.0 , 100.0 ,  5.0 /)
          ZMMD_SOIL(:,5)  = (/ 520.0 ,  75.0 ,  2.5 /)
          ZMMD_SOIL(:,6)  = (/ 520.0 ,  75.0 ,  2.5 /)
          ZMMD_SOIL(:,7)  = (/ 210.0 ,  75.0 ,  2.5 /)
          ZMMD_SOIL(:,8)  = (/ 210.0 ,  50.0 ,  2.5 /)
          ZMMD_SOIL(:,9)  = (/ 125.0 ,  50.0 ,  1.0 /)
          ZMMD_SOIL(:,10) = (/ 100.0 ,  10.0 ,  1.0 /)
          ZMMD_SOIL(:,11) = (/ 100.0 ,  10.0 ,  0.5 /)
          ZMMD_SOIL(:,12) = (/ 100.0 ,  10.0 ,  0.5 /)
          ZSIGMA_SOIL(:,1)  = (/1.6 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,2)  = (/1.6 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,3)  = (/1.6 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,4)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,5)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,6)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,7)  = (/1.7 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,8)  = (/1.7 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,9)  = (/1.7 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,10) = (/1.8 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,11) = (/1.8 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,12) = (/1.8 , 1.8 , 1.8 /)
          ZPCENT_SOIL(:,1)  = (/1.0 , 0.00 , 0.00 /)
          ZPCENT_SOIL(:,2)  = (/0.90 , 0.10 , 0.00 /)
          ZPCENT_SOIL(:,3)  = (/0.80 , 0.20 , 0.00 /)
          ZPCENT_SOIL(:,4)  = (/0.50 , 0.35 , 0.15 /)
          ZPCENT_SOIL(:,5)  = (/0.45 , 0.40 , 0.15 /)
          ZPCENT_SOIL(:,6)  = (/0.35 , 0.50 , 0.15 /)
          ZPCENT_SOIL(:,7)  = (/0.30 , 0.50 , 0.20 /)
          ZPCENT_SOIL(:,8)  = (/0.30 , 0.50 , 0.20 /)
          ZPCENT_SOIL(:,9)  = (/0.20 , 0.50 , 0.30 /)
          ZPCENT_SOIL(:,10) = (/0.65 , 0.00 , 0.35 /)
          ZPCENT_SOIL(:,11) = (/0.60 , 0.00 , 0.40 /)
          ZPCENT_SOIL(:,12) = (/0.50 , 0.00 , 0.50 /)
        ELSEIF (KDDUST==4 .OR. KDDUST==5) THEN
          ZMMD_SOIL(:,1)  = (/1000.0 , 100.0 , 10.0 /)
          ZMMD_SOIL(:,2)  = (/ 690.0 , 100.0 , 10.0 /)
          ZMMD_SOIL(:,3)  = (/ 520.0 , 100.0 ,  5.0 /)
          ZMMD_SOIL(:,4)  = (/ 520.0 , 100.0 ,  5.0 /)
          ZMMD_SOIL(:,5)  = (/ 520.0 ,  75.0 ,  2.5 /)
          ZMMD_SOIL(:,6)  = (/ 520.0 ,  75.0 ,  2.5 /)
          ZMMD_SOIL(:,7)  = (/ 210.0 ,  75.0 ,  2.5 /)
          ZMMD_SOIL(:,8)  = (/ 210.0 ,  50.0 ,  2.5 /)
          ZMMD_SOIL(:,9)  = (/ 125.0 ,  50.0 ,  1.0 /)
          ZMMD_SOIL(:,10) = (/ 100.0 ,  10.0 ,  1.0 /)
          ZMMD_SOIL(:,11) = (/ 100.0 ,  10.0 ,  0.5 /)
          ZMMD_SOIL(:,12) = (/ 100.0 ,  10.0 ,  0.5 /)
          ZSIGMA_SOIL(:,1)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,2)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,3)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,4)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,5)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,6)  = (/1.6 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,7)  = (/1.7 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,8)  = (/1.7 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,9)  = (/1.7 , 1.7 , 1.8 /)
          ZSIGMA_SOIL(:,10) = (/1.8 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,11) = (/1.8 , 1.8 , 1.8 /)
          ZSIGMA_SOIL(:,12) = (/1.8 , 1.8 , 1.8 /)
          ZPCENT_SOIL(:,1)  = (/0.90 , 0.10 , 0.00 /)
          ZPCENT_SOIL(:,2)  = (/0.60 , 0.30 , 0.10 /)
          ZPCENT_SOIL(:,3)  = (/0.60 , 0.30 , 0.10 /)
          ZPCENT_SOIL(:,4)  = (/0.50 , 0.35 , 0.15 /)
          ZPCENT_SOIL(:,5)  = (/0.45 , 0.40 , 0.15 /)
          ZPCENT_SOIL(:,6)  = (/0.35 , 0.50 , 0.15 /)
          ZPCENT_SOIL(:,7)  = (/0.30 , 0.50 , 0.20 /)
          ZPCENT_SOIL(:,8)  = (/0.30 , 0.50 , 0.20 /)
          ZPCENT_SOIL(:,9)  = (/0.20 , 0.50 , 0.30 /)
          ZPCENT_SOIL(:,10) = (/0.65 , 0.00 , 0.35 /)
          ZPCENT_SOIL(:,11) = (/0.60 , 0.00 , 0.40 /)
          ZPCENT_SOIL(:,12) = (/0.50 , 0.00 , 0.50 /)

        ENDIF
        
        ZAEROSIZE_DD = (/1.0E-08 ,2.0E-08 , 4.0E-08 , &
            & 8.0E-08 , 1.6E-07 , 3.2E-07 , &
            & 6.4E-07 , 1.28E-06 ,      &
            & 2.56E-06 , 5.12E-06 , 10.24E-06 , &
            & 20.48E-06 ,  40.96E-06 /)

        ZRHO_DD = YAEROSOLS(1)%RDENSITY

        PFLXDUST(KIDIA:KFDIA,:) = 0.0_JPRB
        ZFSOIL(:,:) = 0.0_JPRB
        
        !-----------------------------------------------------------------------
        
        !*       0.1   SOIL PROPERTIES
        !              --------------------------
        
        IF (KDDUST == 2 ) THEN ! case P. Nabat formulation (CMIP6 version)
        
          !-----------------------------------------------------------------------------------------
          ! -------------- Definition of soil textures
          !
          ! CLAY : CLAY >= 0.40 SILT < 0.40 SAND < 0.45
          ! SANDY CLAY : CLAY >= 0.36 SAND >= 0.45
          ! SILTY CLAY : CLAY >= 0.40 SILT >= 0.40
          ! SILT : SILT >= 0.8 CLAY < 0.12
          ! SAND : SAND >= 0.3*CLAY + 0.87
          ! SANDY CLAY LOAM : CLAY >= 0.28  CLAY < 0.36 SAND >= 0.45 | CLAY >= 0.20 CLAY < 0.28 SILT < 0.28
          ! SILTY CLAY LOAM : CLAY >= 0.28 CLAY < 0.40 SAND < 0.20
          ! CLAY LOAM : CLAY >= 0.28 CLAY < 0.40 SAND >= 0.20 SAND < 0.45
          ! SILT LOAM : SILT >= 0.8 CLAY >= 0.12 | SILT >= 0.5 SILT < 0.8 CLAY < 0.28
          ! LOAMY SAND : SAND >= CLAY + 0.7 SAND < 0.3*CLAY + 0.87
          ! SANDY LOAM : SAND >= 0.52  CLAY < 0.20 | SAND >= (0.5 - CLAY)  CLAY < 0.07
          ! LOAM : CLAY >= 0.20 CLAY < 0.28 SILT >= 0.28 SILT < 0.5 | SAND >= (0.5 - CLAY)  CLAY < 0.20
        
          ZFTEX(:,:) = 0.0_JPRB
        
          DO JL=KIDIA,KFDIA
          !
            ZSAND=PAERDEP(JL)/100
            ZCLAY=PAERSCC(JL)/100
            ZSILT = 1.0_JPRB - ZCLAY - ZSAND
        
            IF (ZSILT <= 0.) THEN
              ZSILT = 0.0_JPRB
            ENDIF
            !
            !tous les cas CLAY >= 0.28 sont traités
            IF ( ZCLAY >= 0.28_JPRB ) THEN
              IF ( ZSAND >= 0.45_JPRB ) THEN
                IF (ZCLAY >= 0.36_JPRB ) THEN     ! Sandy Clay
                  ZFTEX(9,JL)=1.0_JPRB
                ELSE                              ! Sandy Clay Loam
                  ZFTEX(6,JL)=1.0_JPRB
                ENDIF
              ELSEIF ( ZCLAY >= 0.40_JPRB ) THEN
                IF ( ZSILT >= 0.40_JPRB ) THEN    ! Silty Clay
                  ZFTEX(10,JL)=1.0_JPRB
                ELSE                              ! Clay
                  ZFTEX(11,JL)=1.0_JPRB
                ENDIF
              ELSEIF (ZSAND >= 0.20_JPRB ) THEN   ! Clay Loam
                ZFTEX(8,JL)=1.0_JPRB
              ELSE                                ! Silty Clay Loam
                ZFTEX(7,JL)=1.0_JPRB
              ENDIF
            ENDIF
          ! les cas ZSILT >= 0.5 .AND. PCLAY < 0.28 sont traités
          ! les cas ZSILT < 0.5 .AND. PCLAY >= 0.20 sont traités
          ! ne sont pas traités les cas PCLAY < 0.20 .AND. ZSILT < 0.5 => SAND >= 0.3
            IF ( ZSILT >= 0.8_JPRB .AND. ZCLAY < 0.12_JPRB ) THEN ! Silt
              ZFTEX(12,JL)=1.0_JPRB
            ELSEIF ( ZCLAY < 0.28_JPRB ) THEN    ! ( clay est forcément < 0.28 )
              IF ( ZSILT >= 0.5_JPRB ) THEN      ! Silt Loam
                ZFTEX(4,JL)=1.0_JPRB
              ELSEIF ( ZCLAY >= 0.20_JPRB ) THEN
                IF ( ZSILT >= 0.28_JPRB ) THEN   ! Loam
                  ZFTEX(5,JL)=1.0_JPRB
                ELSE                           ! Sandy Clay Loam
                  ZFTEX(6,JL)=1.0_JPRB
                ENDIF
              ENDIF
            ENDIF
          ! les cas SAND >= 0.87 sont traités: silt < 0.13, clay < 0.1 => entrent dans les cas non encore traités
          ! les cas SAND >= 0.7 => CLAY < 0.15 , SILT < 0.3 sont traités ) => ""
          ! les cas SAND >=0.52 .AND. CLAY < 0.20 sont traités SILT < 0.48 => ""
          ! les cas restants considèrement pclay < 0.20, sand >= 0.5 - pclay => sand >=  0.3
          ! => clay + sand >= 0.5 => silt < 0.5
            IF ( ZSAND >= (0.3_JPRB*ZCLAY + 0.87_JPRB) ) THEN   ! Sand
              ZFTEX(1,JL)=1.0_JPRB
            ELSEIF ( ZSAND >= (ZCLAY + 0.7_JPRB) ) THEN    ! Loamy Sand
              ZFTEX(2,JL)=1.0_JPRB
            ELSEIF ( ZSAND >= 0.52_JPRB .AND. ZCLAY < 0.20_JPRB ) THEN ! Sandy Loam
              ZFTEX(3,JL)=1.0_JPRB
            ELSEIF ( ZSAND >= (0.5_JPRB - ZCLAY) ) THEN
              IF ( ZCLAY < 0.07_JPRB ) THEN                   ! Sandy Loam
                ZFTEX(3,JL)=1.0_JPRB
              ELSEIF ( ZCLAY < 0.20_JPRB ) THEN               ! Loam
                ZFTEX(5,JL)=1.0_JPRB
              ENDIF
            ENDIF
          !
          ENDDO
        
        ENDIF  ! KDDUST==2
        
        IF (KDDUST == 4 .OR. KDDUST == 5) THEN ! case P. Nabat formulation (CAMS43 version)
        
          ZFTEX(:,:) = 0.0_JPRB
        
          DO JL=KIDIA,KFDIA
            !
            ZSAND=PAERDEP(JL) / 100._JPRB
            ZCLAY=PAERSCC(JL) / 100._JPRB
            ZSILT = 1.0_JPRB - ZCLAY - ZSAND
        
            IF (ZSILT <= 0._JPRB) THEN
              ZSILT = 0.0_JPRB
            ENDIF
          !
        
            IF (ZSILT + 1.5_JPRB*ZCLAY < 0.15_JPRB) THEN                                                      ! Sand
              ZFTEX(1,JL)=1.0_JPRB
            ELSEIF (ZSILT + 2._JPRB*ZCLAY < 0.30_JPRB) THEN                                                   ! Loamy sand
              ZFTEX(2,JL)=1.0_JPRB
            ELSEIF ((ZCLAY>=0.07_JPRB .AND. ZCLAY<0.20_JPRB .AND. ZSAND>0.52_JPRB .AND. ZSILT+2._JPRB*ZCLAY>=0.30_JPRB) &
              &  .OR. (ZCLAY< 0.07_JPRB .AND. ZSILT<0.50_JPRB .AND. ZSILT+2._JPRB*ZCLAY>=0.30_JPRB)) THEN     ! Sandy Loam
              ZFTEX(3,JL)=1.0_JPRB
            ELSEIF ((ZSILT>=0.50_JPRB .AND. ZCLAY>=0.12_JPRB .AND. ZCLAY<0.27_JPRB) .OR. &
              & (ZSILT>=0.50_JPRB .AND. ZSILT<0.80_JPRB  .AND. ZCLAY<0.12_JPRB)) THEN                         ! Silt Loam
              ZFTEX(4,JL)=1.0_JPRB
            ELSEIF (ZSILT>=0.80_JPRB .AND. ZCLAY<0.12_JPRB) THEN                                              ! Silt
              ZFTEX(5,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.07_JPRB .AND. ZCLAY<0.27_JPRB .AND. ZSILT>=0.28_JPRB &
              & .AND. ZSILT<0.50_JPRB .AND. ZSAND<=0.52_JPRB) THEN ! Loam
              ZFTEX(6,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.20_JPRB .AND. ZCLAY<0.35_JPRB .AND. ZSILT<0.28_JPRB .AND. ZSAND>0.45_JPRB) THEN  ! Sandy Clay Loam
              ZFTEX(7,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.27_JPRB .AND. ZCLAY<0.40_JPRB .AND. ZSAND<=0.20_JPRB) THEN                       ! Silty Clay Loam
              ZFTEX(8,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.27_JPRB .AND. ZCLAY<0.40_JPRB .AND. ZSAND>0.20_JPRB .AND. ZSAND<=0.45_JPRB) THEN ! Clay Loam
              ZFTEX(9,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.35_JPRB .AND. ZSAND>0.45_JPRB) THEN                                              ! Sandy Clay
              ZFTEX(10,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.40_JPRB .AND. ZSILT>=0.40_JPRB) THEN                                             ! Silty Clay
              ZFTEX(11,JL)=1.0_JPRB
            ELSEIF (ZCLAY>=0.40_JPRB .AND. ZSAND<=0.45_JPRB .AND. ZSILT<0.40_JPRB) THEN                       ! Clay
              ZFTEX(12,JL)=1.0_JPRB
            ENDIF
          ENDDO
        
        ENDIF  !KDDUST=4 OR 5
        
        
        
        !-----------------------------------------------------------------------
        
        !*       1.0  DUST EMISSIONS FROM CMIP6
        !              --------------------------
          
        IF (KDDUST == 2 ) THEN 
          
           ! Calculation of ZDP_ARRAY
           ZDP_ARRAY(1)=0.0001_JPRB  !cm
           DO JS=2,INSOIL_DD
             ZDP_ARRAY(JS) = ZDP_ARRAY(JS-1)*EXP(0.0460517018598807_JPRB)
           ENDDO
          
           ! Emission distribution (Kok, 2011)
           ZTOTV = 0.0_JPRB
           DO JN=1,ISIZE_DD
             ZRWI = (ZAEROSIZE_DD(JN)+ZAEROSIZE_DD(JN+1))/2.0*1.E6
             ZFRAC(JN) = ZRWI/ZCV_KOK*(1+ERF(LOG(ZRWI/ZDS_KOK)/SQRT(2.0)/LOG(ZSIGMAS_KOK)))*EXP(-(ZRWI/ZLAMBDA_KOK)**3) !Kok (2011)
             ZTOTV = ZTOTV + ZFRAC(JN)
           ENDDO
          
           DO JN=1,ISIZE_DD
             ZFRAC(JN) = ZFRAC(JN)/ZTOTV
             if ( ZFRAC(JN) < 1.E-9 ) ZFRAC(JN) = 0.0_JPRB
           ENDDO
          
           !-- PRELIMINARY FOR TESTING
          
           !-- surface source of dust is assumed if fraction of bare soil is > 0.1 and/or UVis albedo > 0.11
           !                                        10m wind > threshold = f(soil wetness, particle radius)
           !                                        no snow, and below the wilting point (0.171)
           DO JL=KIDIA,KFDIA
             LLDUST(JL)=.TRUE.
          
             !-------------------------------------------------------------------
             ! ---  MB95 - Modifications by P.Nabat (10/2013)
          
             !      IF (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND. PBARESOIL(JL) >= 0.90_JPRB &
             !      Sensibilite test adopted
             IF (LDUST_SFXSOIL) THEN
               LLDDACT(JL)= (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND. PBARESOIL(JL) >= 0.50_JPRB &
               &   .AND. ZSWI1(JL) <= 0.5_JPRB )
             ELSE
               LLDDACT(JL)= (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND. PBARESOIL(JL) >= 0.50_JPRB &
               &   .AND. PWS1(JL) <= 2._JPRB * ZRWPWP )
                 !.AND. PHSDFOR(JL) <= 50._JPRB
             ENDIF
          
             ! ZFEFF = correction factor for the effect of surface roughness
             IF (LLDDACT(JL)) THEN
               ZFEFF(JL)=1._JPRB - (LOG(MIN(0.5E-4_JPRB,PZ0M(JL))/ZZ0SAL)/LOG(0.35_JPRB*(0.001_JPRB/ZZ0SAL)**0.8_JPRB))
          
               ! ZFW = factor for the effect of soil moisture
               IF (LDUST_SFXSOIL) THEN
                 ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB, ZSWI1(JL) ) )
               ELSE
                 ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB, (PWS1(JL)-ZRWPWP)/(ZRWSAT-ZRWPWP) ) )
               ENDIF
               ZSWETN2 = 0.0015_JPRB*PAERSCC(JL)**2 + 0.17_JPRB*PAERSCC(JL)
               IF (ZSWETN <= ZSWETN2) THEN
                 ZFW(JL) = 1._JPRB
               ELSE
                 ZFW(JL) = (1._JPRB + 1.21_JPRB*(ZSWETN-ZSWETN2)**0.68_JPRB)**0.5_JPRB   ! Fecan et al., 1999
               ENDIF
             ENDIF
           ENDDO
          
           ! ZUSTARTS = ideal minimum threshold friction velocity
           DO JS = 1,INSOIL_DD
             ZRE = 1331.647_JPRB*ZDP_ARRAY(JS)**1.561228+0.38194
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 ZKK = SQRT(2*YDCST_ACLIB%RG*ZRHO_DD*100*ZDP_ARRAY(JS)/PRHOS(JL) &
                      & *(1+0.006/(ZRHO_DD*YDCST_ACLIB%RG*(2*100*ZDP_ARRAY(JS))**2.5)))
                 IF (ZRE > 10) THEN
                   ZUSTARTS(JL,JS)=0.129*ZKK*(1-0.0858*exp(-0.0617*(ZRE-10)))
                 ELSE
                   ZUSTARTS(JL,JS)=0.129*ZKK/(1.928*(ZRE**0.092)-1)**0.5
                 ENDIF
                 ZUSTART(JL,JS) = ZUSTARTS(JL,JS)/ZFEFF(JL)*ZFW(JL)
               ENDIF
             ENDDO
           ENDDO
          
           ! Determination of zsrel
           DO JNT = 1 , INATS_DD
             ZSS(:) = 0.0_JPRB
             ZSTOTAL = 0.0_JPRB
             DO JS = 1 , INSOIL_DD          !soil size segregatoin no
               DO JNM = 1 , INSOILMODE_DD       !soil mode = 3
                 IF ( (ZPCENT_SOIL(JNM,JNT) > ZEPSDD) .AND. (ZSIGMA_SOIL(JNM,JNT)/=0.0) ) then
                   ZXK = ZPCENT_SOIL(JNM,JNT)/(SQRT(2*YDCST_ACLIB%RPI)*LOG(ZSIGMA_SOIL(JNM,JNT)))
                   ZXL = ((LOG(ZDP_ARRAY(JS))-LOG(ZMMD_SOIL(JNM,JNT)*1.E-4))**2)/(2.0*(LOG(ZSIGMA_SOIL(JNM,JNT)))**2)
                   ZXM = ZXK*EXP(-ZXL)
                 ELSE
                   ZXM = 0.0_JPRB
                 ENDIF
                 ZXN = ZRHO_DD*(2.0_JPRB/3.0_JPRB)*(ZDP_ARRAY(JS)/2.0_JPRB)
                 ZDELDP = 0.0460517018598807_JPRB        ! zdp_array(2)-zdp_array(1) zss(nsoil)
                 ZSS(JS) = ZSS(JS) + (ZXM*ZDELDP/ZXN)
               ENDDO
               ZSTOTAL = ZSTOTAL + ZSS(JS)
             ENDDO
             DO JS = 1 , INSOIL_DD
               IF (ZSTOTAL > 0.0_JPRB) THEN
                 ZSREL(JNT,JS) = ZSS(JS)/ZSTOTAL
               ENDIF
             ENDDO
           ENDDO                    ! soil types
          
           ! Saltation flux (Marticorena and Bergametti, 1995)
           DO JL=KIDIA,KFDIA
             IF (LLDDACT(JL)) THEN
               ZUTMIN = (15._JPRB/(100._JPRB*ZRVONKAR*ZFEFF(JL)))*LOG(10._JPRB/PZ0M(JL))
               IF ( PWIND(JL) >= ZUTMIN )  THEN
                 ZUSTAR(JL) = 100._JPRB*PWIND(JL) + 0.3_JPRB*(PWIND(JL)-ZUTMIN)**2
               ELSE
                 ZUSTAR(JL) = 100._JPRB*PWIND(JL)
               ENDIF
             ENDIF
           ENDDO
          
           DO JS=1,INSOIL_DD
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 IF (ZUSTAR(JL) /= 0._JPRB) THEN
                   ZUTH = ZUSTART(JL,JS)/ZUSTAR(JL)
                 ELSE
                   ZUTH = 0._JPRB
                 ENDIF
                 ZFDP1 = ZUSTAR(JL)**3 * MAX(0.0_JPRB,(1.0_JPRB-ZUTH**2))
                 ZFDP2 = (1._JPRB+ZUTH)*(1.0E-5)*PRHOS(JL)/YDCST_ACLIB%RG
                 IF (PAERSCC(JL) <= 17.0_JPRB) THEN  !before correction in Eastern Africa: 20.0
                   ZALPHAPROP = YPAER_ACLIB%RCOEFDUST*1.0E-6*10**(0.134*PAERSCC(JL)-6.0)
                 ELSE
                   !ZALPHAPROP = YPAER_ACLIB%RCOEFDUST*1.0E-6*10**(-0.1*PAERSCC(JL)-1.2)
                   ZALPHAPROP = YPAER_ACLIB%RCOEFDUST*1.0E-6*10**(-0.09*PAERSCC(JL)-2.19)
                 ENDIF
          
                 DO JNT=1,INATS_DD
                   ZFSOIL(JNT,JL) = ZFSOIL(JNT,JL) + ZALPHAPROP*ZSREL(JNT,JS)*ZFDP1*ZFDP2
                 ENDDO
               ENDIF
             ENDDO
           ENDDO
          
           ! Calculate vertical fluxes
           ZRSFROWT(:,KIDIA:KFDIA,:) = 0.0_JPRB
           DO JAER=1,KAERO
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 DO JN=1,ISIZE_DD
                   DO JNT=1,INATS_DD
                     ZRSFROWSUB(JNT,JN,JL) = ZFSOIL(JNT,JL)*ZFRAC(JN)
                     ZRWI = (ZAEROSIZE_DD(JN)+ZAEROSIZE_DD(JN+1))/2.0*1.E6
                     IF (ZRWI >= ZTRSIZE_DD(JAER) .AND. ZRWI < ZTRSIZE_DD(JAER+1)) THEN
                       ZRSFROWT(JNT,JL,JAER) = ZRSFROWT(JNT,JL,JAER) + ZRSFROWSUB(JNT,JN,JL)
                     ENDIF
                   ENDDO
                 ENDDO
               ENDIF
             ENDDO ! JL
          
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 DO JNT=1,INATS_DD
                   PFLXDUST(JL,JAER) = PFLXDUST(JL,JAER) + ZRSFROWT(JNT,JL,JAER) * ZFTEX(JNT,JL) &
                                     &  * PBARESOIL(JL) * (1._JPRB-PSNS(JL)) * PLSM(JL)
                 ENDDO
               ENDIF
             ENDDO
          
           ENDDO ! enddo KAERO
          
        ENDIF ! case KDDUST == 2
          
        !-----------------------------------------------------------------------
        
        !*       2.0  UPDATED DUST EMISSIONS (developed in CAMS43)
        !              --------------------------
        
        IF (KDDUST == 4 .OR. KDDUST == 5 ) THEN 
          
           ! Calculation of ZDP_ARRAY
           ZDELDP = 0.0460517018598807_JPRB        ! zdp_array(2)-zdp_array(1) zss(nsoil)dln(Dp)
           ZDP_ARRAY(1)=0.0001_JPRB  ! in cm
           DO JS=2,INSOIL_DD
                ZDP_ARRAY(JS) = ZDP_ARRAY(JS-1)*exp(ZDELDP)
           ENDDO
          
           ! Determination of ZSREL=dSrel in MB95
           ! equations 29-30-31 in MB95
           ! ZSS=dS(Dp)
           DO JNT = 1 , INATS_DD
             ZSS(:) = 0.0_JPRB
             ZSTOTAL = 0.0_JPRB
             DO JS = 1 , INSOIL_DD          !soil size segregation
               DO JNM = 1 , INSOILMODE_DD       !soil mode = 3
                 IF ( ZPCENT_SOIL(JNM,JNT)>0 .AND. ZSIGMA_SOIL(JNM,JNT)>0 ) THEN
                   ZXK = ZPCENT_SOIL(JNM,JNT)/(sqrt(2*YDCST_ACLIB%RPI)*LOG(ZSIGMA_SOIL(JNM,JNT)))
                   ZXL =((LOG(ZDP_ARRAY(JS))-LOG(ZMMD_SOIL(JNM,JNT)*1.E-4))**2)/(2.0*(LOG(ZSIGMA_SOIL(JNM,JNT)))**2)
                   ZXM = ZXK*EXP(-ZXL)
                 ELSE
                   ZXM = 0.0_JPRB
                 ENDIF
                 ZXN = 2.0*YDCST_ACLIB%RPI*(ZDP_ARRAY(JS)*1.E4)**2
                 ZSS(JS) = ZSS(JS) + (ZXM*ZDELDP*ZXN)
               ENDDO
               ZDDP=((EXP(ZDELDP))**(JS-1))*(EXP(ZDELDP)-1)*ZDP_ARRAY(1)
               ZSTOTAL = ZSTOTAL + ZSS(JS)*ZDDP
             ENDDO
             DO JS = 1 , INSOIL_DD
               IF (ZSTOTAL > 0.0 ) THEN
                 ZSREL(JNT,JS) = ZSS(JS)/ZSTOTAL
               ENDIF
             ENDDO
           ENDDO                    ! soil types   
          
           ! Emission distribution (Kok, 2011)
           ZTOTV = 0.0_JPRB
           DO JN=1,ISIZE_DD
             ZRWI = (ZAEROSIZE_DD(JN)+ZAEROSIZE_DD(JN+1))/2.0*1.E6
             ZFRAC(JN) = ZRWI/ZCV_KOK*(1+ERF(LOG(ZRWI/ZDS_KOK)/SQRT(2.0)/LOG(ZSIGMAS_KOK)))*EXP(-(ZRWI/ZLAMBDA_KOK)**3)
             ZTOTV = ZTOTV + ZFRAC(JN)
           ENDDO
          
           DO JN=1,ISIZE_DD
             ZFRAC(JN) = ZFRAC(JN)/ZTOTV
           ENDDO
          
           !-- PRELIMINARY FOR TESTING
          
           !-- surface source of dust is assumed if fraction of bare soil is > 0.1
           !   and/or UVis albedo > 0.11, 10m wind > threshold = f(soil wetness, particle radius)
           !   no snow, and below the wilting point (0.171)
        
           DO JL=KIDIA,KFDIA
        
              LLDUST(JL)=.TRUE.
        
              !-------------------------------------------------------------------
              ! ---  MB95 - Modifications by P.Nabat (10/2013)
          
              IF (LDUST_SFXSOIL) THEN
                  LLDDACT(JL)= (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND. PBARESOIL(JL) >= 0.70_JPRB &
                   & .AND. ZSWI1(JL) <= 0.5_JPRB .AND. PWIND(JL) > 0._JPRB)
              ELSE
                  LLDDACT(JL)= (PLSM(JL) >= 0.99_JPRB .AND. PSNS(JL) == 0.0_JPRB .AND. PBARESOIL(JL) >= 0.70_JPRB &
                   & .AND. PWS1(JL) <= 2._JPRB * ZRWPWP  .AND. PWIND(JL) > 0._JPRB)
                    !.AND. PHSDFOR(JL) <= 50._JPRB
              ENDIF
           ENDDO
          
           IF (KDDUST == 5) THEN
              DO JS=1,INSOIL_DD
                ! Calculation of Z0s : roughness length of smooth erodible surfaces
                ! Equ (18) in Marticorena et al. (1997)
                ZZ0S(JS) = (ZDP_ARRAY(JS)*1.E-2)/30.
                DO JL=KIDIA,KFDIA
                  IF (LLDDACT(JL)) THEN
                    ! ZFEFFS = correction factor for the effect of surface roughness
                    ! Equ (3) in Zakey et al. (2006)
                    ! Equ (20) in Marticorena and Bergametti (1995)
                    ! Equ (2) in Laurent et al. (2008)
                    ! Equ (4) in Menut et al. (2013) with error in units (Z0s should be in m)
                    !ZFEFFS(JL,JS)=MIN(1.0,MAX(0.001,1-(LOG(PZ0M(JL)/(200.*ZZ0S(JS)))/LOG(0.35*(0.1/ZZ0S(JS))**0.8))))
                    ZFEFFS(JL,JS) = MIN(1.0_JPRB, MAX(0.001_JPRB, 1._JPRB - &
                           & (LOG(PZ0M(JL)/(100._JPRB*ZZ0S(JS)))/LOG(0.35_JPRB*(0.1_JPRB/ZZ0S(JS))**0.8_JPRB))))
                  ENDIF
                ENDDO
              ENDDO
           ELSE   ! KDDUST==4
              DO JL=KIDIA,KFDIA
                ! ZFEFF = correction factor for the effect of surface roughness
                IF (LLDDACT(JL)) THEN
                  ZFEFF(JL) = MIN(1.0_JPRB, MAX(0.001_JPRB,1._JPRB - &
                           & (LOG(0.0002_JPRB/ZZ0SAL)/LOG(0.35_JPRB*(0.1_JPRB/ZZ0SAL)**0.8_JPRB))))
                ENDIF
              ENDDO
           ENDIF  ! KDDUST 4/5
        
           DO JL=KIDIA,KFDIA
              IF (LLDDACT(JL)) THEN
                ! ZFW = factor for the effect of soil moisture = Wetu*T/DryU*T
                ! ZSWETN=w "total soil moisture" in Fecan et al. (1999)
                ! w="prognostic surface volumetric soil moisture calculated by BATS in Zakey et al.(2006)
                IF (LDUST_SFXSOIL) THEN
                  ZSWETN = MAX(0.001_JPRB, 100*(ZSWI1(JL)*(ZRWSAT-ZRWPWP)+ZRWPWP) )
                ELSE
                  ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB, 100*PWS1(JL) ) )
                ENDIF
                ! ZSWETN2=w' in Fecan et al. (1999) : empirical equation (14/15)
                ! w'=residual soil moisture (Laurent et al. 2008)
                ZSWETN2 = 0.0014_JPRB*PAERSCC(JL)**2+0.17_JPRB*PAERSCC(JL)   ! 0.0015 in Zakey et al. (2006) ?
                IF (ZSWETN <= ZSWETN2) THEN
                  ZFW(JL)=1._JPRB
                ELSE
                  ZFW(JL)=(1._JPRB+1.21_JPRB*(ZSWETN-ZSWETN2)**0.68_JPRB)**0.5_JPRB   ! Fecan et al.(1999) equ (15)
                ENDIF
              ENDIF
           ENDDO
          
           ! Calculation of ZUSTARTS = ideal minimum threshold friction velocity
           !   and ZUSTART = threshold wind friction velocity
           DO JS = 1,INSOIL_DD
             ! Marticorena and Bergametti (1995) :
             ! ZRE= Reynolds number (equation 5)
             ZRE = 1331.647*ZDP_ARRAY(JS)**1.561228+0.38194    ! With ZDP_ARRAY in cm
          
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 ! Marticorena and Bergametti (1995), ZKK=K in equation (4)
                 ZKK = (ZRHO_DD/PRHOS(JL) * YDCST_ACLIB%RG * ZDP_ARRAY(JS)*1.E-2_JPRB * & ! RRHO_DD and PRHOS in kg.m-3
                     & (1._JPRB + 0.006_JPRB/(0.001_JPRB * ZRHO_DD * 100._JPRB*YDCST_ACLIB%RG * & 
                     & (ZDP_ARRAY(JS))**2.5)))**0.5  !Put RG in cm.s-2 ???
                 IF (ZRE > 10._JPRB) THEN
                   ! MB95 equation (7)
                   ZUSTARTS(JL,JS)=0.129_JPRB*ZKK*(1._JPRB-0.0858_JPRB*EXP(-0.0617_JPRB*(ZRE-10._JPRB)))
                 ELSE
                   ! MB95 equation (6)
                   ZUSTARTS(JL,JS)=0.129_JPRB*ZKK/((1.928_JPRB*(ZRE**0.092_JPRB)-1._JPRB)**0.5_JPRB)
                 ENDIF
                 ! MB95 equation (10/17/21)
                 ! Zakey et al. (2006) equ (1)
                 ! Laurent et al. (2008) equ (1)
                 ! Menut et al. (2013) equ (3)
                 IF (KDDUST == 4) THEN
                   ZUSTART(JL,JS) = ZUSTARTS(JL,JS)/ZFEFF(JL)*ZFW(JL)
                 ELSE
                   ZUSTART(JL,JS) = ZUSTARTS(JL,JS)/ZFEFFS(JL,JS)*ZFW(JL)
                 ENDIF
               ENDIF
             ENDDO
           ENDDO
          
           DO JL=KIDIA,KFDIA
            IF (LLDDACT(JL)) THEN
                ZUSTAR(JL) = PWIND(JL)*ZRVONKAR/LOG(10._JPRB/PZ0M(JL))
            ENDIF
           ENDDO
          
           DO JS=1,INSOIL_DD
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL) .AND. ZUSTAR(JL)>ZUSTART(JL,JS)) THEN
                 ZUTH = ZUSTART(JL,JS)/ZUSTAR(JL)
                 ZFDP1 = ZUSTAR(JL)**3._JPRB * MAX(0._JPRB,(1.0_JPRB-ZUTH**2._JPRB))
                 ZFDP2 = (1._JPRB+ZUTH)*PRHOS(JL)/YDCST_ACLIB%RG
          
                 ! Conversion from horizontal to vertical flux (sandblasting)
                 ! Gilette (1979) - equ (47) in MB95       
                 IF (PAERSCC(JL) <= 20.0_JPRB) THEN  !before correction in Eastern Africa: 20.0
                   ZALPHAPROP = YPAER_ACLIB%RCOEFDUST*10._JPRB**(0.134_JPRB*PAERSCC(JL)-6.0_JPRB)
                 ELSEIF (PAERSCC(JL) <= 45.0_JPRB) THEN
                   ZALPHAPROP = YPAER_ACLIB%RCOEFDUST*10._JPRB**(-0.067_JPRB*PAERSCC(JL)-3.0_JPRB)
                 ELSE
                   ZALPHAPROP = 10._JPRB**(-0.018_JPRB*PAERSCC(JL)-5.18_JPRB)
                 ENDIF
          
                 IF (PAERDEP(JL) > 0._JPRB) THEN
                   DO JNT=1,INATS_DD
                     ZFSOIL(JNT,JL) = ZFSOIL(JNT,JL) + ZALPHAPROP * ZSREL(JNT,JS) * ZFDP1 * ZFDP2 * &
                                    & ((EXP(ZDELDP))**(JS-1))*(EXP(ZDELDP)-1._JPRB)*ZDP_ARRAY(1)
                   ENDDO
                 ENDIF
               ENDIF
             ENDDO
           ENDDO
        
           ! Calculate vertical fluxes (Kok distribution)
           ZRSFROWT(:,:,:) = 0.0_JPRB
        
           DO JAER=1,KAERO
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 DO JN=1,ISIZE_DD
                   DO JNT=1,INATS_DD
                     ZRSFROWSUB(JNT,JN,JL) = ZFSOIL(JNT,JL)*ZFRAC(JN)
                     ZRWI = (ZAEROSIZE_DD(JN)+ZAEROSIZE_DD(JN+1))/2.0_JPRB*1.E6_JPRB
                     IF (ZRWI >= ZTRSIZE_DD(JAER) .AND. ZRWI < ZTRSIZE_DD(JAER+1)) THEN
                       ZRSFROWT(JNT,JL,JAER) = ZRSFROWT(JNT,JL,JAER) + ZRSFROWSUB(JNT,JN,JL)
                     ENDIF
                   ENDDO
                 ENDDO
               ENDIF
             ENDDO
          
             DO JL=KIDIA,KFDIA
               IF (LLDDACT(JL)) THEN
                 DO JNT=1,INATS_DD
                   PFLXDUST(JL,JAER) = PFLXDUST(JL,JAER) + &
                     & ZRSFROWT(JNT,JL,JAER)*ZFTEX(JNT,JL)*PBARESOIL(JL)*(1._JPRB-PSNS(JL))*PLSM(JL)
                 ENDDO
               ENDIF
             ENDDO
           ENDDO ! enddo NTYPAER(2)
        ENDIF ! case KDDUST == 2 or 4
        
        !-----------------------------------------------------------------------
        
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_EMIS_DUST:TACTIC_EMIS_DUST',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
    END SUBROUTINE TACTIC_EMIS_DUST

END MODULE M_TACTIC_EMIS_DUST
