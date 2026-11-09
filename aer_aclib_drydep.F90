MODULE M_AER_ACLIB_DRYDEP
 
    CONTAINS

    SUBROUTINE AER_ACLIB_DRYDEP(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                          & CDAERO_CODE, KAERO, &
                          & KAER_START, KIDIA, KFDIA, &
                          & KLUOUT, &
                          & KKA, KKU, KKL, & 
                          & KLON, KLAT, KLEV, PTSTEP, &
                          & PAP, PAPH, PAPHI, PDELP, PRHODREF, &
                          & PCI, PLSM, PSNS, &
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
        ! PAP (KLON,KLEV)         : FULL-LEVEL PRESSURE           (Pa)
        ! PAPH (KLON,0:KLEV)      : HALF-LEVEL PRESSURE           (Pa)
        ! PAPHI(KLON,0:KLEV)      : GEOPOTENTIAL ON HALF-LEVELS 
        ! PDELP (KLON,KLEV)       : DELTA PRESSURE                (Pa)
        ! PRHODEF(KLON,KLEV)      : AIR VOLUMIC MASS
        ! PCI(KLON)               : FRACTION OF SEA-ICE
        ! PLSM(KLON)              : LAND-SEA MASK (note: treats lakes as sea)
        ! PSNS(KLON)              : SNOW DEPTH

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
        USE M_TACTIC_DRYDEP,      ONLY : TACTIC_DRYDEP

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
        REAL(KIND=JPRB),    INTENT(IN)  :: PAP(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PAPH(KLON,0:KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PAPHI(KLON,0:KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PDELP(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PRHODREF(KLON,KLEV)  
        REAL(KIND=JPRB),    INTENT(IN)  :: PCI(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PLSM(KLON) 
        REAL(KIND=JPRB),    INTENT(IN)  :: PSNS(KLON)

       ! INPUT/OUTPUT VARIABLES 
        REAL(KIND=JPRB),     INTENT(INOUT) :: PAER(KLON,KLEV,KAERO)   
        REAL(KIND=JPRB),     INTENT(INOUT) :: PTENAER(KLON,KLEV,KAERO)   
        REAL(KIND=JPRB),     INTENT(INOUT) :: PSFAER(KLON,KAERO)   
        TYPE(TAEROUT_ACLIB), INTENT(INOUT) :: YAEROUT_ACLIB 

       ! LOCAL VARIABLES
        ! General
        !INTEGER(KIND=JPIM) :: JAER, JL
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        ! For tactic
        !REAL(KIND=JPRB)    :: ZFAERO(KLON, KAERO), ZVDEP(KLON)


!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB_DRYDEP:AER_ACLIB_DRYDEP',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
 
        SELECT CASE (CDAERO_CODE)

            CASE ("MOCAGE")
                ! To do

            CASE ("tactic")

                !! Dry deposition !!
                !ZFAERO(KIDIA:KFDIA,1:KAERO) = 0.0_JPRB

                ! WRITE(KLUOUT,*) "AER_ACLIB_DRYDEP : TACTIC - ",YAEROSOLS(1)%CAER_NAME
                CALL TACTIC_DRYDEP (YDCST_ACLIB, YAEROSOLS, &
                    &  KIDIA , KFDIA , KLON , KLEV, KAERO, &
                    &  KAER_START, KLUOUT, PTSTEP, &
                    &  PAER, PTENAER, &
                    &  PAPH, PAP, & 
                    &  PDELP, PAPHI, PRHODREF(:,KLEV) , &
                    &  PLSM, PCI, PSNS, &
                    &  PSFAER, YAEROUT_ACLIB)
                    !&  ZFAERO, ZVDEP)

                ! No need => PSFAER and YAEROUT_ACLIB directly modified in TACTIC_DRYDEP
                !DO JAER=1,KAERO
                !  WRITE(KLUOUT,*) "AER_ACLIB_DRYDEP : TACTIC - ",JAER," - ",YAEROSOLS(JAER)%CAER_NAME
                !  DO JL=KIDIA,KFDIA
                !    YAEROUT_ACLIB%PDRYDEP(JL,JAER+KAER_START-1) = ZFAERO(JL,JAER)
                !    PSFAER(JL,JAER) = PSFAER(JL,JAER) + ZFAERO(JL,JAER)
                !  ENDDO
                !ENDDO


            CASE ("MNHDUST")
                ! To do
  
            CASE DEFAULT
                ! WRITE(KLUOUT,*) "aer_aclib_drydep : ",TRIM(CDAERO_CODE), "not yet implemented"
                !STOP

        END SELECT
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB_DRYDEP:AER_ACLIB_DRYDEP',1,ZHOOK_HANDLE)    
!------------------------------------------------------------------------------

    END SUBROUTINE AER_ACLIB_DRYDEP

END MODULE M_AER_ACLIB_DRYDEP
