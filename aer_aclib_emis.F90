MODULE M_AER_ACLIB_EMIS
 
    CONTAINS

    SUBROUTINE AER_ACLIB_EMIS(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                          & CDAERO_CODE, KAERO, &
                          & KAER_START, KIDIA, KFDIA, &
                          & KLUOUT, &
                          & KKA, KKU, KKL, & 
                          & KLON, KLAT, KLEV, PTSTEP, &
                          & PRH, PRHODREF, &
                          & PWIND, PSST, PWS1, PCI, PBARESOIL, PLSM, PCLK, PSNS, &
                          & PAERDEP, PAERSCC, PZ0M, &
                          & PAER, PTENAER, PSFAER, &
                          & YAEROUT_ACLIB) 

        ! INPUT:
        ! -------
        ! YDCST_ACLIB             : STRUCTURE FOR CONSTANTS
        ! YPAER_ACLIB             : STRUCTURE FOR GENERAL PARAMETERS OF AEROSOLS
        ! YAEROSOLS               : STRUCTURE FOR PARAMETERS OF EACH AEROSOL
        ! CDAERO_CODE             : AEROSOL CODE TO BE USED FOR SURFACE FLUXES
        ! KAERO                   : NUMBER OF AEROSOLS FOR THE CURRENT PROCESS
        ! KAER_START              : FIRST INDEX OF AEROSOL BIN/MODE TO BE USED FOR THIS PROCESS
        ! KIDIA , KFDIA           : BEGIN/END of LON BLOCK
        ! KLUOUT                  : UNIT OF OUTPUT FILE
        ! KKA                      : Near ground array index (is an unphysical level if JPVEXT!=0)
        ! KKU                      : Uppest atmosphere array index (is an unphysical level if JPVEXT!=0)
        ! KKL                     : Order of the vertical levels:
        !                         1: as for Méso-NH, levels are numbered from ground to space
        !                        -1: as for AROME/MOCAGE/ARPIFS, levels are numbered from space to ground
        ! KLON , KLAT             : NUMBER OF LONGITUDES / LATITUDES
        ! KLEV                    : NUMBER OF VERTICAL LEVELS
        ! PTSTEP                  : TIME STEP 
        ! PRH(KLON,KLEV)          : RELATIVE HUMIDITY (0-1)
        ! PRHODEF(KLON,KLEV)      : AIR VOLUMIC MASS
        ! PWIND                   : Wind speed (average of horizontal wind speed) (m/s)
        ! PSST                    : Sea surface temperature (K)
        ! PWS1                    : soil humidity
        ! PCI(KLON)               : FRACTION OF SEA-ICE
        ! PBARESOIL(KLON)         : BARE SOIL FRACTION
        ! PLSM(KLON)              : LAND-SEA MASK (note: treats lakes as sea)
        ! PCLK(KLON)              : LAKE COVER
        ! PSNS(KLON)              : SNOW DEPTH
        ! PAERDEP(KLON)           : SAND CONTENT (%)
        ! PAERSCC(KLON)           : CLAY CONTENT (%)
        ! PZ0M(KLON)              : ROUGHNESS LENGTH (m)

        ! INPUT/OUTPUT:
        ! -------
        ! PAER(KLON,KLEV,KAERO)   : AEROSOL CONCENTRATION
        ! PTENAER(KLON,KLEV,KAERO): AEROSOL TENDENCY
        ! PSFAER(KLON,KAERO)      : SURFACE AEROSOL FLUX

        ! OUTPUT:
        ! -------
        ! YAEROUT_ACLIB           : STRUCTURE FOR AEROSOL OUTPUT

        USE PARKIND1,             ONLY : JPIM, JPRB
        USE YOMHOOK,              ONLY : LHOOK,   DR_HOOK
        USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
        USE MODD_PARAM_AER_ACLIB, ONLY : TPARAM_AER_ACLIB, TAEROSOL, TAEROUT_ACLIB
        USE M_TACTIC_EMIS_SALT,   ONLY : TACTIC_EMIS_SALT
        USE M_TACTIC_EMIS_DUST,   ONLY : TACTIC_EMIS_DUST

        IMPLICIT NONE

       ! STRUCTURES
        TYPE(TCST_ACLIB),                INTENT(IN)  :: YDCST_ACLIB
        TYPE(TPARAM_AER_ACLIB),          INTENT(IN)  :: YPAER_ACLIB
        TYPE(TAEROSOL),DIMENSION(KAERO), INTENT(IN)  :: YAEROSOLS
 
       ! INPUT VARIABLES
        CHARACTER(LEN=*),   INTENT(IN)  :: CDAERO_CODE
        INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA , KFDIA 
        INTEGER(KIND=JPIM), INTENT(IN)  :: KAERO
        INTEGER(KIND=JPIM), INTENT(IN)  :: KAER_START
        INTEGER(KIND=JPIM), INTENT(IN)  :: KLON , KLAT, KLEV
        INTEGER(KIND=JPIM), INTENT(IN)  :: KLUOUT,KKA, KKU, KKL
        REAL(KIND=JPRB),    INTENT(IN)  :: PTSTEP
        REAL(KIND=JPRB),    INTENT(IN)  :: PRH(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PRHODREF(KLON,KLEV)  
        REAL(KIND=JPRB),    INTENT(IN)  :: PWIND(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PSST(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PWS1(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PCI(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PBARESOIL(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PLSM(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PCLK(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PSNS(KLON)
        REAL(KIND=JPRB),    INTENT(IN)  :: PAERDEP(KLON)
        REAL(KIND=JPRB),    INTENT(IN)  :: PAERSCC(KLON)
        REAL(KIND=JPRB),    INTENT(IN)  :: PZ0M(KLON)

       ! INPUT/OUTPUT VARIABLES 
        REAL(KIND=JPRB),     INTENT(INOUT) :: PAER(KLON,KLEV,KAERO)   
        REAL(KIND=JPRB),     INTENT(INOUT) :: PTENAER(KLON,KLEV,KAERO)   
        REAL(KIND=JPRB),     INTENT(INOUT) :: PSFAER(KLON,KAERO)   
        TYPE(TAEROUT_ACLIB), INTENT(INOUT) :: YAEROUT_ACLIB 

       ! LOCAL VARIABLES
        ! General
        INTEGER(KIND=JPIM) :: JAER, JL
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        ! For tactic
        REAL(KIND=JPRB)    :: ZFAERO(KLON, KAERO)
        INTEGER(KIND=JPIM) :: IDUST
        LOGICAL            :: LLDUST_SFXSOIL, LLSST


!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB_EMIS:AER_ACLIB_EMIS',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
 
        SELECT CASE (CDAERO_CODE)

            CASE ("MOCAGE")
                ! To do

            CASE ("tactic")

                !! Emissions (only for dust and sea-salt, other aerosols are emitted in SURFEX)

                IF ( YAEROSOLS(1)%CAERFAM_NAME == "SALT" .AND. &
                    & YAEROSOLS(KAERO)%CAERFAM_NAME == "SALT" .AND. &
                    & YAEROSOLS(1)%NBINMOD == KAERO) THEN
                    ! Sea-salt
                    SELECT CASE (YAEROSOLS(1)%CEMIS)
                        CASE ("grythe-sst")
                            LLSST = .TRUE.
                        CASE ("grythe")
                            LLSST = .FALSE.
                        CASE DEFAULT
                            WRITE(KLUOUT,*) "Unknown emission scheme for dust aerosols"
                    END SELECT 
                    ! WRITE(KLUOUT,*) "AER_ACLIB_EMIS TACTIC - LLSST ",LLSST
                    CALL TACTIC_EMIS_SALT(YDCST_ACLIB, YAEROSOLS, &
                        & KIDIA, KFDIA, KLON, KAERO, KLUOUT, LLSST,&
                        & PCI, PLSM, PCLK, PWIND, PSST, PRH(:,KLEV), ZFAERO)

                ELSEIF ( YAEROSOLS(1)%CAERFAM_NAME == "DUST" .AND. &
                    & YAEROSOLS(KAERO)%CAERFAM_NAME == "DUST" .AND. &
                    & YAEROSOLS(1)%NBINMOD == KAERO) THEN
                    ! Dust             
                    ! WRITE(KLUOUT,*) "ACLIB TACTIC EMIS - AEROSOL",YAEROSOLS(1)%CAER_NAME
                    SELECT CASE (YAEROSOLS(1)%CEMIS)
                        CASE ("nddust2")
                            IDUST = 2
                            LLDUST_SFXSOIL = .FALSE.
                        CASE ("nddust4")
                            IDUST = 4
                            LLDUST_SFXSOIL = .FALSE.
                        CASE ("nddust5")
                            IDUST = 5
                            LLDUST_SFXSOIL = .FALSE.
                        CASE ("nddust2-sfx")
                            IDUST = 2
                            LLDUST_SFXSOIL = .TRUE.
                        CASE ("nddust4-sfx")
                            IDUST = 4
                            LLDUST_SFXSOIL = .TRUE.
                        CASE ("nddust5-sfx")
                            IDUST = 5
                            LLDUST_SFXSOIL = .TRUE.
                        CASE DEFAULT
                            WRITE(KLUOUT,*) "Unknown emission scheme for dust aerosols"
                    END SELECT 
                    !WRITE(KLUOUT,*) "AER_ACLIB_EMIS TACTIC - IDUST ",IDUST
                    CALL TACTIC_EMIS_DUST(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                        & KIDIA, KFDIA, KLON, KAERO, KLUOUT, &
                        & IDUST, LLDUST_SFXSOIL, &
                        & PBARESOIL, PLSM, PSNS, PAERDEP, PAERSCC, PZ0M, PWIND, PRHODREF(:,KLEV), PWS1, &
                        & ZFAERO)

                ENDIF

                DO JAER=1,KAERO
                !   WRITE(KLUOUT,*) "AER_ACLIB_EMIS : TACTIC - ",JAER," - ",YAEROSOLS(JAER)%CAER_NAME, "Mean : ", SUM(ZFAERO(:,JAER)) / SIZE(ZFAERO, DIM=1)
                  DO JL=KIDIA,KFDIA
                    YAEROUT_ACLIB%PEMIS(JL,JAER+KAER_START-1) = ZFAERO(JL,JAER)
                    PSFAER(JL,JAER) = PSFAER(JL,JAER) - ZFAERO(JL,JAER)
                  ENDDO
                ENDDO


            CASE ("MNHDUST")
                ! To do
  
            CASE DEFAULT

                WRITE(KLUOUT,*) "aer_aclib_emis : ",TRIM(CDAERO_CODE), "not yet implemented"
                !STOP

        END SELECT
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB_EMIS:AER_ACLIB_EMIS',1,ZHOOK_HANDLE)    
!------------------------------------------------------------------------------

    END SUBROUTINE AER_ACLIB_EMIS

END MODULE M_AER_ACLIB_EMIS
