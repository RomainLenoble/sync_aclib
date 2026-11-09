MODULE M_AER_ACLIB_SEDIM
 
    CONTAINS

    SUBROUTINE AER_ACLIB_SEDIM(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
                          & CDAERO_CODE, KAERO, KMODE,&
                          & KAER_START, KIDIA, KFDIA, &
                          & KLUOUT, &
                          & KKA, KKU, KKL, & 
                          & KLON, KLAT, KLEV, PTSTEP, &
                          & PAP, PAPH, PDELP, PDZZ, PTP, PVERVEL, PRH, PRHODREF, &
                          & CDMNH_PROGRAM, &
                          & PAER, PTENAER, & 
                          & YAEROUT_ACLIB) 

        ! INPUT:
        ! -------
        ! YDCST_ACLIB             : STRUCTURE FOR CONSTANTS
        ! YPAER_ACLIB             : STRUCTURE FOR GENERAL PARAMETERS OF AEROSOLS
        ! YAEROSOLS               : STRUCTURE FOR PARAMETERS OF EACH AEROSOL
        ! CDAERO_CODE             : AEROSOL CODE TO BE USED FOR SEDIMENTATION
        ! CDMNH_PROGRAM           : CALLING PROGRAM FOR MESO-NH
        ! KAERO                   : NUMBER OF AEROSOLS FOR THE CURRENT PROCESS
        ! KAER_START              : FIRST INDEX OF AEROSOL BIN/MODE TO BE USED FOR THIS PROCESS
        ! KIDIA , KFDIA           : BEGIN/END of LON BLOCK
        ! KLUOUT                  : UNIT OF OUTPUT FILE
        ! KKA                     : Near ground array index (is an unphysical level if JPVEXT!=0)
        ! KKU                     : Uppest atmosphere array index (is an unphysical level if JPVEXT!=0)
        ! KKL                     : Order of the vertical levels:
        !                         1: as for Méso-NH, levels are numbered from ground to space
        !                        -1: as for AROME/MOCAGE/ARPIFS, levels are numbered from space to ground
        ! KLON , KLAT             : NUMBER OF LONGITUDES / LATITUDES
        ! KLEV                    : NUMBER OF VERTICAL LEVELS
        ! PTSTEP                  : TIME STEP 
        ! PAP (KLON,KLEV)         : FULL-LEVEL PRESSURE           (Pa)
        ! PAPH (KLON,0:KLEV)      : HALF-LEVEL PRESSURE           (Pa)
        ! PDELP (KLON,KLEV)       : DELTA PRESSURE                (Pa)
        ! PDZZ(KLON,KLEV)         : HEIGHT OF LAYERS              (m)
        ! PTP(KLON,KLEV)          : FULL-LEVEL TEMPERATURE        (K)
        ! PVERVEL(KLON,KLEV)      : FULL-LEVEL VERTICAL VELOCITY  (Pa s-1)
        ! PRH(KLON,KLEV)          : RELATIVE HUMIDITY             (0-1)
        ! PRHODEF(KLON,KLEV)      : AIR VOLUMIC MASS

        ! INPUT/OUTPUT:
        ! -------
        ! PAER(KLON,KLEV,KAERO)   : AEROSOL CONCENTRATION
        ! PTENAER(KLON,KLEV,KAERO): AEROSOL TENDENCY

        ! OUTPUT:
        ! -------
        ! YAEROUT_ACLIB           : STRUCTURE FOR AEROSOL OUTPUT

        USE PARKIND1,             ONLY : JPIM, JPRB
        USE YOMHOOK,              ONLY : LHOOK, DR_HOOK
        USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
        USE MODD_PARAM_AER_ACLIB, ONLY : TPARAM_AER_ACLIB, TAEROSOL, TAEROUT_ACLIB
        USE M_MF_MOCAGE_SEDIM,    ONLY : MF_MOCAGE_SEDIM
        USE M_TACTIC_SEDIMNT,     ONLY : TACTIC_SEDIMNT
        USE MODI_SEDIM_DSLT,      ONLY : SEDIM_DSLT
        USE MODI_CH_AER_SEDIM_n,  ONLY:  CH_AER_SEDIM_n 

        IMPLICIT NONE

       ! STRUCTURES
        TYPE(TCST_ACLIB),                INTENT(IN)  :: YDCST_ACLIB
        TYPE(TPARAM_AER_ACLIB),          INTENT(IN)  :: YPAER_ACLIB
        TYPE(TAEROSOL),DIMENSION(KAERO), INTENT(IN)  :: YAEROSOLS
 
       ! INPUT VARIABLES
        CHARACTER(LEN=*),   INTENT(IN)  :: CDAERO_CODE
        INTEGER(KIND=JPIM), INTENT(IN)  :: KAERO
        INTEGER(KIND=JPIM), INTENT(IN)  :: KAER_START
        INTEGER(KIND=JPIM), INTENT(IN)  :: KMODE
        INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA , KFDIA 
        INTEGER(KIND=JPIM), INTENT(IN)  :: KLUOUT, KKA, KKU, KKL
        INTEGER(KIND=JPIM), INTENT(IN)  :: KLON , KLAT, KLEV
        REAL(KIND=JPRB),    INTENT(IN)  :: PTSTEP
        REAL(KIND=JPRB),    INTENT(IN)  :: PAP(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PAPH(KLON,0:KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PDELP(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PDZZ(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PTP(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PVERVEL(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PRH(KLON,KLEV)   
        REAL(KIND=JPRB),    INTENT(IN)  :: PRHODREF(KLON,KLEV)  
        CHARACTER(LEN=*),   INTENT(IN)  :: CDMNH_PROGRAM

       ! INPUT/OUTPUT VARIABLES 
        REAL(KIND=JPRB),     INTENT(INOUT) :: PAER(KLON,KLEV,KAERO)   
        REAL(KIND=JPRB),     INTENT(INOUT) :: PTENAER(KLON,KLEV,KAERO)
        TYPE(TAEROUT_ACLIB), INTENT(INOUT) :: YAEROUT_ACLIB 

       ! LOCAL VARIABLES
        ! General
        INTEGER(KIND=JPIM) :: JAER, JK, JL
        REAL(KIND=JPRB)  :: ZDUMMY(KLON)           
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        ! For mocage        
        REAL(KIND=JPRB),DIMENSION(KAERO)   :: ZDIAMED,ZDENSAER
        CHARACTER(LEN=15),DIMENSION(KAERO) :: CLAERNAM
        ! For tactic
        REAL(KIND=JPRB)  :: ZTAERA(KLON,KLEV), ZTAERZ(KLON,KLEV), ZSEDIMV(KLON,KLEV)   
        REAL(KIND=JPRB)  :: ZFAER(KLON)
        REAL(KIND=JPRB)  :: ZFAER4(KLON,KLEV,KAERO)
        ! For mesonh
        INTEGER(KIND=JPIM)  :: ISPCARBSOA     ! NSP+NCARB+NSOA
        REAL(KIND=JPRB)  :: ZTHT(KLON,KLEV)   ! Dry potential temperature

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB_SEDIM:AER_ACLIB_SEDIM',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
 
        SELECT CASE (CDAERO_CODE)

            CASE ("MOCAGE")
             
                ! Gather aerosol propeties for call sedim
                DO JAER=1,KAERO
                    ZDIAMED(JAER) = YAEROSOLS(JAER)%RDIAMED
                    ZDENSAER(JAER) = YAEROSOLS(JAER)%RDENSITY
                    CLAERNAM(JAER) = YAEROSOLS(JAER)%CAER_NAME
                    DO JK=1,KLEV
                      DO JL=KIDIA,KFDIA
                        ZFAER4(JL,JK,JAER)=PAER(JL,JK,JAER)                 
                      ENDDO
                    ENDDO
                ENDDO
 
                CALL MF_MOCAGE_SEDIM (YDCST_ACLIB, &
                    &       KIDIA,KFDIA,1,1,1,KLON,1,1, &
                    &       KLEV,KAERO,PTSTEP, &
                    &       PAP,PDELP,PTP,ZFAER4,ZDIAMED,ZDENSAER, &
                    &       PRH,CLAERNAM, KLUOUT)

                IF ( YPAER_ACLIB%LWORK_IN_TEN ) THEN

                    DO JAER=1,KAERO
                      DO JK=1,KLEV
                        DO JL=KIDIA,KFDIA
                          PTENAER(JL,JK,JAER)=(ZFAER4(JL,JK,JAER)-PAER(JL,JK,JAER))/PTSTEP      
                        ENDDO
                      ENDDO
                    ENDDO

                ELSE

                    DO JAER=1,KAERO
                      DO JK=1,KLEV
                        DO JL=KIDIA,KFDIA
                          PAER(JL,JK,JAER)=ZFAER4(JL,JK,JAER)                 
                        ENDDO
                      ENDDO
                    ENDDO
                 
                ENDIF

            CASE ("tactic")

                DO JAER=1,KAERO
                      
                    IF (YAEROSOLS(JAER)%CSEDIM == "tactic") THEN

                        !WRITE(KLUOUT,*) "aclib_sedim : YAEROSOLS(",JAER,")%RSEDIMVEL=",YAEROSOLS(JAER)%RSEDIMVEL
                        DO JK=1,KLEV
                          DO JL=KIDIA,KFDIA
                            ZTAERA(JL,JK)=PTENAER(JL,JK,JAER)
                            ZTAERZ(JL,JK)=0.
                            ZSEDIMV(JL,JK)=YAEROSOLS(JAER)%RSEDIMVEL
                          ENDDO
                        ENDDO

                        !WRITE(KLUOUT,*) "aclib_sedim : before PTENAER(1,KLEV,",JAER,")=",PTENAER(1,KLEV,JAER)
                        CALL TACTIC_SEDIMNT (YDCST_ACLIB, &
                            &  KIDIA , KFDIA , KLON , KLEV, ZSEDIMV, &
                            &  ZDUMMY   , ZDUMMY  , &
                            &  PTSTEP, PTP   , PAP, PAPH, PVERVEL, PRHODREF, & 
                            &  PAER(:,:,JAER)  , ZTAERA, &
                            &  ZTAERZ, ZFAER ) 

                        ! Output to be saved
                        DO JL=KIDIA,KFDIA
                            YAEROUT_ACLIB%PSEDIM(JL,JAER+KAER_START-1) = ZFAER(JL) 
                        ENDDO

                        IF ( YPAER_ACLIB%LWORK_IN_TEN ) THEN
                            DO JK=1,KLEV
                              DO JL=KIDIA,KFDIA
                                PTENAER(JL,JK,JAER)=ZTAERZ(JL,JK)
                              ENDDO
                            ENDDO
                        ELSE
                            DO JK=1,KLEV
                              DO JL=KIDIA,KFDIA
                                PAER(JL,JK,JAER)=MAX(PAER(JL,JK,JAER)+ZTAERZ(JL,JK)*PTSTEP,YDCST_ACLIB%REPSI)
                              ENDDO
                            ENDDO
                        ENDIF
                         
                        !WRITE(KLUOUT,*) "aclib_sedim : after PTENAER(1,KLEV,",JAER,")=",PTENAER(1,KLEV,JAER)
 
                    ENDIF

                ENDDO

            CASE ("MNHDSLT")

                ZTHT(:,:) = PTP(:,:) / (PAP(:,:) / YDCST_ACLIB%RP00)**(YDCST_ACLIB%RRD/YDCST_ACLIB%RCPD)
                CALL SEDIM_DSLT(YDCST_ACLIB, YAEROSOLS, KKA, KKU, KKL, KLON, KLAT, KLEV, KAERO, KMODE, &
                                ZTHT, PTSTEP, PRHODREF,  PAP, PDZZ, PAER )

            CASE ("MNHORILAM")

                 ZTHT(:,:) = PTP(:,:) / (PAP(:,:) / YDCST_ACLIB%RP00)**(YDCST_ACLIB%RRD/YDCST_ACLIB%RCPD)
                 ISPCARBSOA = (YPAER_ACLIB%NANTH_F1 - YPAER_ACLIB%NANTH_F0 + 1) / KMODE
                 CALL CH_AER_SEDIM_n(YDCST_ACLIB,YAEROSOLS, CDMNH_PROGRAM,&
                                     PTSTEP, KKA,KKU,KKL,KLON,KLAT,KLEV,KAERO,KMODE, ISPCARBSOA, &
                                     PAER, ZTHT, PRHODREF, PAP, PDZZ, YAEROUT_ACLIB%PSEDIM)

            CASE DEFAULT
                WRITE(KLUOUT,*) "aer_aclib_sedim : ",TRIM(CDAERO_CODE), "not yet implemented"
                !STOP

        END SELECT
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_AER_ACLIB_SEDIM:AER_ACLIB_SEDIM',1,ZHOOK_HANDLE)    
!------------------------------------------------------------------------------

    END SUBROUTINE AER_ACLIB_SEDIM

END MODULE M_AER_ACLIB_SEDIM
