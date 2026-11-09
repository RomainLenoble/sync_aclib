MODULE M_AER_ACLIB
 
    CONTAINS

    SUBROUTINE AER_ACLIB(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                    & CDAERO_PROCESS, CDAERO_CODE, &
                    & KAERO, KMODE, KCHEM,  &
                    & KAER_START, KIDIA, KFDIA, &
                    & KLUOUT, KVERB, &
                    & KKA, KKU, KKL, & 
                    & KLON, KLAT, KLEV, PTSTEP, &
                    & PAER, PTENAER, &
                    & YAEROUT_ACLIB, &
                    & YPCHEM_ACLIB, YCHEMS, &
                    & PTIME, &
                    & PAP, PAPH, PDELP, PRHODREF, PTP, &
                    & PVERVEL, PRH, PDZZ, &
                    & PAPHI, PSST, PWS1, PCI, PBARESOIL, PLSM, &
                    & PCLK, PSNS, PAERDEP, PAERSCC, PZ0M, PWIND, &
                    & PEVAP3D, PVMASSMIN, PRGAER, &
                    & PLCVPFLX, PLSTPFLX, PSCVPFLX, PSSTPFLX, PSURF, &
                    & PCLWATS, PCLCOV, PCLWATL, PGELAT, PRV, PRC, &
                    & CDMNH_PROGRAM, PCCT, PCRT, PSEA, PTOWN, PAERDENSITY3D, &
                    & LDCALLRAD, KTOP, &
                    & PSFAER, PCHEM, PTENCHEM, YCHORILAM_PACKED)

        ! INPUTS:
        ! -------
        ! CDAERO_PROCESS          : AEROSOL PROCESSUS
        ! CDAERO_CODE             : AEROSOL CODE TO BE USED FOR THIS PROCESS
        ! KAERO                   : NUMBER OF AEROSOLS FOR THE CURRENT PROCESS
        ! KMODE                   : NUMBER OF MODES FOR THE CURRENT PROCESS
        ! KCHEM                   : NUMBER OF CHEMICAL SPECIES
        ! KAER_START              : FIRST INDEX OF AEROSOL BIN/MODE TO BE USED FOR THIS PROCESS
        ! KIDIA , KFDIA           : BEGIN/END of LON BLOCK
        ! KLUOUT                  : UNIT OF OUTPUT FILE
        ! KVERB                   : VERBOSITY LEVEL
        ! KKA                     : Near ground array index (is an unphysical level if JPVEXT!=0)
        ! KKU                     : Uppest atmosphere array index (is an unphysical level if JPVEXT!=0)
        ! KKL                     : Order of the vertical levels:
        !                         1: as for Méso-NH, levels are numbered from ground to space
        !                        -1: as for AROME/MOCAGE/ARPIFS, levels are numbered from space to ground
        ! KLON , KLAT             : NUMBER OF LONGITUDES / LATITUDES
        ! KLEV                    : NUMBER OF VERTICAL LEVELS
        ! PTSTEP                  : TIME STEP 

        ! INPUTS OPTIONAL:
        ! -------
        ! CDMNH_PROGRAM           : CALLING PROGRAM FOR MESO-NH
        ! PTIME                   : SECONDS SINCE 00 UTC
        ! PAP (KLON,KLEV)         : FULL-LEVEL PRESSURE           (Pa)
        ! PAPH (KLON,0:KLEV)      : HALF-LEVEL PRESSURE           (Pa)
        ! PAPHI(KLON,0:KLEV)      : GEOPOTENTIAL ON HALF-LEVELS 
        ! PDELP (KLON,KLEV)       : DELTA PRESSURE                (Pa)
        ! PRHODREF(KLON,KLEV)     : AIR VOLUMIC MASS
        ! PTP(KLON,KLEV)          : FULL-LEVEL TEMPERATURE        (K)
        ! PVERVEL(KLON,KLEV)      : FULL-LEVEL VERTICAL VELOCITY  (Pa s-1)
        ! PRH(KLON,KLEV)          : RELATIVE HUMIDITY             (0-1)
        ! 
        ! CASE SEDIM
        ! -------
        ! PDZZ(KLON,KLEV)         : HEIGHT OF LAYERS              (m)
        !
        ! CASE EMIS AND DRYDEP
        ! -------
        ! PAPHI(KLON,0:KLEV)      : GEOPOTENTIAL ON HALF-LEVELS 
        ! PSST                    : Sea surface temperature (K)
        ! PWS1                    : Soil humidity
        ! PCI(KLON)               : FRACTION OF SEA-ICE
        ! PBARESOIL(KLON)         : BARE SOIL FRACTION
        ! PLSM(KLON)              : LAND-SEA MASK (note: treats lakes as sea)
        ! PCLK(KLON)              : LAKE COVER
        ! PTOWN(KLON)             : TOWN FRACTION
        ! PSNS(KLON)              : SNOW DEPTH
        ! PAERDEP(KLON)           : SAND CONTENT (%)
        ! PAERSCC(KLON)           : CLAY CONTENT (%)
        ! PZ0M(KLON)              : ROUGHNESS LENGTH (m)
        ! PWIND(KLON)             : SURFACE WIND SPEED           (m/s)   
        !
        ! CASE DRYDEP
        ! -------
        !
        ! CASE SCAV
        ! -------
        ! PEVAP3D(KLON,KLEV)      : Instantaneous 3D Rain Evaporation flux (KG/KG/S)
        ! PVMASSMIN               : Minimum mass to transfer between dry mass or in-cloud droplets
        ! PRGAER                  : Number median diameter [um}
        ! PAERDENSITY3D           : Density in 3D for each mode (??)
        ! PLCVPFLX(KLON,0:KLEV)   : LIQUID CONVECTIVE PRECIPITATION FLUX (kg.m-2.s-1)
        ! PLSTPFLX(KLON,0:KLEV)   : LIQUID STRATIFORM PRECIPITATION FLUX (kg.m-2.s-1)
        ! PSCVPFLX(KLON,0:KLEV)   : SOLID  CONVECTIVE PRECIPITATION FLUX (kg.m-2.s-1)
        ! PSSTPFLX(KLON,0:KLEV)   : SOLID  STRATIFORM PRECIPITATION FLUX (kg.m-2.s-1)
        ! PSURF                   : SURFACE OF THE MESH           (m2) 
        ! PCLWATS(KLON,KLEV)      : SOLID  WATER CONTENT         (kg/kg)
        !
        ! CASES SCAV & TRANSFO
        ! -------
        ! PCLCOV(KLON,KLEV)       : CLOUD COVER                   (0-1)
        ! PCLWATL(KLON,KLEV)      : LIQUID WATER CONTENT         (kg/kg)
        ! PCCT(KLON,KLEV)         : CLOUD LIQUID WATER CONCENTRATION  (m-3)
        ! PCRT(KLON,KLEV)         : RAIN WATER CONCENTRATION          (m-3)
        ! PRV(KLON,KLEV)          : WATER VAPOR MIXING RATIO     (kg/kg)
        ! PRC(KLON,KLEV)          : LIQUID WATER MIXING RATIO    (kg/kg)
        !
        ! CASE TRANSFO
        ! -------
        ! PGELAT                  : LATITUDE
        !
        ! CASE DIAGS
        ! -------
        ! LDCALLRAD                : IS RADIATION CALLED AT THE PRESENT TIME STEP
        ! KTOP (KLON)            : LEVEL OF CLOUD TOP

        ! INPUT/OUTPUT:
        ! -------
        ! PAER(KLON,KLEV,KAERO)    : AEROSOL CONCENTRATION
        ! PTENAER(KLON,KLEV,KAERO) : AEROSOL TENDENCY
        ! YAEROUT_ACLIB            : STRUCTURE FOR AEROSOL OUTPUT
        ! YCHORILAM_PACKED         : STRUCTURE FOR ORILAM (PACKED) INOUTPUT
        
        ! INPUT/OUTPUT OPTIONAL:
        ! -------
        ! SURFLUX
        ! -------
        ! PSFAER(KLON,KAERO)           : AEROSOL SURFACE FLUX
        !
        ! TRANSFO
        ! -------

        USE PARKIND1,              ONLY : JPIM     ,JPRB
        USE YOMHOOK,               ONLY : LHOOK,   DR_HOOK
        USE M_AER_ACLIB_EMIS,      ONLY : AER_ACLIB_EMIS
        USE M_AER_ACLIB_DRYDEP,    ONLY : AER_ACLIB_DRYDEP
        USE M_AER_ACLIB_SEDIM,     ONLY : AER_ACLIB_SEDIM
        USE M_AER_ACLIB_SCAV,      ONLY : AER_ACLIB_SCAV
        USE M_AER_ACLIB_TRANSFO,   ONLY : AER_ACLIB_TRANSFO
        USE M_AER_ACLIB_DIAGS,     ONLY : AER_ACLIB_DIAGS
        USE MODD_CST_ACLIB,        ONLY : TCST_ACLIB
        USE MODD_PARAM_AER_ACLIB,  ONLY : TPARAM_AER_ACLIB, TAEROSOL, TAEROUT_ACLIB
        USE MODD_PARAM_CHEM_ACLIB, ONLY : TPARAM_CHEM_ACLIB, TCHEM
        USE MODD_CH_ORILAM,        ONLY : TCH_ORILAM

        IMPLICIT NONE

        TYPE(TCST_ACLIB),                INTENT(IN)    :: YDCST_ACLIB
        TYPE(TPARAM_AER_ACLIB),          INTENT(IN)    :: YPAER_ACLIB
        TYPE(TAEROSOL),DIMENSION(KAERO), INTENT(IN)    :: YAEROSOLS

        CHARACTER(LEN=*),                INTENT(IN)    :: CDAERO_PROCESS 
        CHARACTER(LEN=*),                INTENT(IN)    :: CDAERO_CODE
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KIDIA , KFDIA 
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KAERO, KAER_START
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KMODE
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KCHEM
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KLUOUT, KVERB
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KLON , KLAT, KLEV
        INTEGER(KIND=JPIM),              INTENT(IN)    :: KKA, KKU, KKL
        REAL(KIND=JPRB),                 INTENT(IN)    :: PTSTEP
        REAL(KIND=JPRB),                 INTENT(INOUT) :: PAER(KLON,KLEV,KAERO)   
        REAL(KIND=JPRB),                 INTENT(INOUT) :: PTENAER(KLON,KLEV,KAERO) 
    
        TYPE(TAEROUT_ACLIB),             INTENT(INOUT) :: YAEROUT_ACLIB

        TYPE(TPARAM_CHEM_ACLIB),     OPTIONAL,INTENT(IN) :: YPCHEM_ACLIB
        TYPE(TCHEM),DIMENSION(KCHEM),OPTIONAL,INTENT(IN) :: YCHEMS

        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PTIME
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PAP(KLON,KLEV)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PAPH(KLON,0:KLEV)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PDELP(KLON,KLEV)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PRHODREF(KLON,KLEV) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PTP(KLON,KLEV)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PVERVEL(KLON,KLEV)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PRH(KLON,KLEV)  
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PDZZ(KLON,KLEV)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PAPHI(KLON,0:KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PSST(KLON) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PWS1(KLON) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCI(KLON) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PBARESOIL(KLON) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PLSM(KLON) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCLK(KLON) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PSNS(KLON)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PAERDEP(KLON)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PAERSCC(KLON)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PZ0M(KLON)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PWIND(KLON)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PEVAP3D(KLON,KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PVMASSMIN(KLON,KLEV,KMODE)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PAERDENSITY3D(KLON,KLEV,KMODE)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PRGAER(KLON,KLEV,KMODE)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PLCVPFLX(KLON,0:KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PLSTPFLX(KLON,0:KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PSCVPFLX(KLON,0:KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PSSTPFLX(KLON,0:KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PSURF(KLON)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCLWATS(KLON,KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCLCOV(KLON,KLEV) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCLWATL(KLON,KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PRV(KLON,KLEV) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PRC(KLON,KLEV) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PGELAT(KLON) 
        LOGICAL,         OPTIONAL,   INTENT(IN)  :: LDCALLRAD
        INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN)  :: KTOP(KLON)
        CHARACTER(LEN=*),OPTIONAL,   INTENT(IN)  :: CDMNH_PROGRAM 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCCT(KLON,KLEV) 
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PCRT(KLON,KLEV)
        REAL(KIND=JPRB), OPTIONAL,   INTENT(IN)  :: PTOWN(KLON),PSEA(KLON)

        REAL(KIND=JPRB), OPTIONAL,   INTENT(INOUT) :: PSFAER(KLON,KAERO)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(INOUT) :: PCHEM(KLON,KLEV,KCHEM)   
        REAL(KIND=JPRB), OPTIONAL,   INTENT(INOUT) :: PTENCHEM(KLON,KLEV,KCHEM) 
        TYPE(TCH_ORILAM),OPTIONAL,   INTENT(INOUT) :: YCHORILAM_PACKED
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB:AER_ACLIB',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------

        SELECT CASE (CDAERO_PROCESS)
 
            !***** SURFACE FLUXES *******
            CASE ("EMIS")       
 !               WRITE(KLUOUT,*) "AER_ACLIB : CALL EMIS STARTING FROM",YAEROSOLS(1)%CAER_NAME
                CALL AER_ACLIB_EMIS(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, & 
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

            CASE ("DRYDEP")       
!               WRITE(KLUOUT,*) "AER_ACLIB : CALL DRYDEP STARTING FROM",YAEROSOLS(1)%CAER_NAME
                CALL AER_ACLIB_DRYDEP(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, & 
                    & CDAERO_CODE, KAERO, &
                    & KAER_START, KIDIA, KFDIA, &
                    & KLUOUT, &
                    & KKA, KKU, KKL, & 
                    & KLON, KLAT, KLEV, PTSTEP, &
                    & PAP, PAPH, PAPHI, PDELP, PRHODREF, &
                    & PCI, PLSM, PSNS, &
                    & PAER, PTENAER, PSFAER, &
                    & YAEROUT_ACLIB) 

            !***** SEDIMENTATION *******
            CASE ("SEDIM")       
!                WRITE(KLUOUT,*) "AER_ACLIB : CALL SEDIM STARTING FROM",YAEROSOLS(1)%CAER_NAME
                CALL AER_ACLIB_SEDIM(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, & 
                    & CDAERO_CODE, KAERO, KMODE, &
                    & KAER_START, KIDIA, KFDIA, &
                    & KLUOUT, &
                    & KKA, KKU, KKL, & 
                    & KLON, KLAT, KLEV, PTSTEP, &
                    & PAP, PAPH, PDELP, PDZZ, PTP, PVERVEL, PRH, PRHODREF, &
                    & CDMNH_PROGRAM, &
                    & PAER, PTENAER, & 
                    & YAEROUT_ACLIB) 

            !***** SCAVENGING *******
            CASE ("SCAV") 
!                WRITE(KLUOUT,*) "AER_ACLIB : CALL SCAV STARTING FROM",YAEROSOLS(1)%CAER_NAME
                CALL AER_ACLIB_SCAV(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                    & CDAERO_CODE, KAERO, KMODE, &
                    & KAER_START, KIDIA, KFDIA, &
                    & KLUOUT, KKA, KKU, KKL, &
                    & KLON, KLAT, KLEV, PTSTEP, &
                    & PAP, PAPH, PDELP, PDZZ, PTP, PVERVEL, PRHODREF, &
                    & PEVAP3D, PVMASSMIN, PRGAER, &
                    & PLCVPFLX, PLSTPFLX, PSCVPFLX, PSSTPFLX, PSURF, &
                    & PCLCOV, PCLWATL, PCLWATS, &
                    & PCCT, PCRT, PSEA,PTOWN,   &
                    & CDMNH_PROGRAM, PAERDENSITY3D, & 
                    & PAER, PTENAER, YAEROUT_ACLIB) 

            !***** TRANSFORMATION *******
            CASE ("TRANSFO") 
!                WRITE(KLUOUT,*) "AER_ACLIB : CALL TRANSFO STARTING FROM",YAEROSOLS(1)%CAER_NAME
                CALL AER_ACLIB_TRANSFO(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                    & CDAERO_CODE, KAERO, KMODE, &
                    & YPCHEM_ACLIB, YCHEMS, &
                    & KCHEM, &
                    & KIDIA, KFDIA, &
                    & KLUOUT, KVERB, &
                    & KLON,  KLEV, PTSTEP, PTIME, PGELAT, &
                    & PAP, & 
                    & PDELP, &
                    & PTP, &
                    & PRH, PRHODREF, PCLCOV, PCLWATL, PRV, PRC, &
                    & YCHORILAM_PACKED,&
                    & PAER,  PTENAER, &
                    & PCHEM ,  PTENCHEM , &
                    & YAEROUT_ACLIB) 

            !***** DIAGS *******
            CASE ("DIAGS") 
!                WRITE(KLUOUT,*) "AER_ACLIB : CALL DIAGS STARTING FROM",YAEROSOLS(1)%CAER_NAME
                CALL AER_ACLIB_DIAGS(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                    & CDAERO_CODE, KAERO,  &
                    & KIDIA, KFDIA, &
                    & KLUOUT, &
                    & KKA, KKU, KKL, & 
                    & KLON, KLAT, KLEV, PTSTEP, &
                    & LDCALLRAD, KTOP, &
                    & PAP, PAPH, PDELP, PTP,  &
                    & PRH, PRHODREF, &
                    & PAER, PTENAER,  &
                    & YAEROUT_ACLIB) 

            CASE DEFAULT
                WRITE(KLUOUT,*) "AER_ACLIB : UNKNOWN AEROSOL PROCESS - STOP "

        END SELECT

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB:AER_ACLIB',1,ZHOOK_HANDLE)    
!------------------------------------------------------------------------------

    END SUBROUTINE AER_ACLIB

END MODULE M_AER_ACLIB
