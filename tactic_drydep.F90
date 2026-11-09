MODULE M_TACTIC_DRYDEP

    CONTAINS
    
    SUBROUTINE TACTIC_DRYDEP &
        !---inputs
         & ( YDCST_ACLIB, YAEROSOLS, &
         &   KIDIA, KFDIA, KLON, KLEV, KAERO, &
         &   KAER_START, KLUOUT, PTSPHY, &
         &   PAERI, PTAERI, &
         &   PAPH, PAP, &
         &   PDP, PGEOH , PRHO, & 
         &   PLSM, PCI, PSNS, &
        !-- outputs
        ! &   PFAERO , PVDEP)
         &   PSFAER, YAEROUT_ACLIB)
        
        
        !**** TACTIC_DRYDEP* -  ROUTINE FOR PARAMETRIZATION OF DRY DEPOSITION
        
        !      from aer_drydep : Jean-Jacques Morcrette 
        !      following O.Boucher's formulation for LMD-Z
        
        !      Dry deposition is (simply) represented by a modification of the 
        !      instantaneous surface flux by what comes down from layer just 
        !      above the surface
        
        !**   INTERFACE.
        !     ----------
        !          *TACTIC_DRYDEP* IS CALLED FROM *AER_ACLIB_DRYDEP*.
        
        ! INPUTS:
        ! -------
        ! YDCST_ACLIB             : STRUCTURE FOR CONSTANTS
        ! YAEROSOLS               : STRUCTURE FOR PARAMETERS OF EACH AEROSOL
        ! PTAERI(KLON,KLEV,KAERO) : INPUT TENDENCIES          (xx kg s-1)
        
        ! OUTPUTS:
        ! --------
        ! PFAERO(KLON,KAERO)      : SURFACE FLUX              (xx m-2)
        
        !     EXTERNALS.
        !     ----------
        !          NONE
        
        !     MODIFICATIONS.
        !     -------------
        !          JJMorcrette 20110725 maximum deposition speed
        !          SRémy       20160309 SO2 dry deposition velocity from SUMO (same as
        !          CHEM)
        !          SRémy       20160830 dry deposition velocities (except SO2) computed
        !          following Zhang et al 2001.
        !          SRémy       20170117 land classes from JFlemming's branch
        !          SRémy       20170512  externalize diurnal cycle
        
        !     SWITCHES.
        !     --------
        
        !     MODEL PARAMETERS
        !     ----------------
        
        !-----------------------------------------------------------------------
        USE PARKIND1,      ONLY : JPIM, JPRB
        USE YOMHOOK,       ONLY : LHOOK, DR_HOOK
        
        USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
        USE MODD_PARAM_AER_ACLIB, ONLY : TAEROSOL, TAEROUT_ACLIB

        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        
        !*       0.1  ARGUMENTS
        !             ---------
        
        !---input fields
        ! STRUCTURES
        TYPE(TCST_ACLIB),                INTENT(IN)  :: YDCST_ACLIB
        TYPE(TAEROSOL),DIMENSION(KAERO), INTENT(IN)  :: YAEROSOLS

        INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KAER_START
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLUOUT
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 

        REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERI(KLON,KLEV,KAERO)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAERI(KLON,KLEV,KAERO)
   
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,0:KLEV)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV), PGEOH(KLON,0:KLEV)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHO(KLON)
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON), PCI(KLON) 
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNS(KLON) ! Snow depth
        
        
        !---output fields
        !REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAERO(KLON,KAERO)
        !REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDEP(KLON,KAERO)
        REAL(KIND=JPRB),     INTENT(INOUT) :: PSFAER(KLON,KAERO)   
        TYPE(TAEROUT_ACLIB), INTENT(INOUT) :: YAEROUT_ACLIB 
        
        !--- local variables
        INTEGER(KIND=JPIM) :: JL, JAER
        REAL(KIND=JPRB)    :: ZAERO, ZFROC, ZOCEA, ZSICE, ZLAND
        REAL(KIND=JPRB)    :: ZVDEP(KAERO)
        REAL(KIND=JPRB)    :: ZMAXVDRY(KLON), ZDZ(KLON)
        
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_DRYDEP:TACTIC_DRYDEP',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
        !- PFAERO    in unit of xx m-2 s-1 (surface flux)
        !- N.B. Surface emission fluxes are negative upward, so contribution of 
        !  dry deposition is to make the surface flux less negative 
        
        ! PFAERO(KIDIA:KFDIA,1:KAERO) = 0.0_JPRB

        ZVDEP(:)=0._JPRB
        
        DO JL=KIDIA,KFDIA
          ZDZ(JL)=(PAPH(JL,KLEV)-PAPH(JL,KLEV-1)) / (YDCST_ACLIB%RG*PRHO(JL))  
          ZMAXVDRY(JL)= ZDZ(JL) / PTSPHY
        ENDDO
     
        
        DO JAER=1,KAERO
        
            IF (TRIM(YAEROSOLS(JAER)%CDRY) == "1") THEN
                ! Skip species if no deposition velocity configured
              IF (.NOT.(YAEROSOLS(JAER)%RDRYDEPVELOCE == 0._JPRB .AND. &
                      & YAEROSOLS(JAER)%RDRYDEPVELLND == 0._JPRB)) THEN
        
                ! WRITE(KLUOUT,*) "TACTIC_DRYDEP - AER",JAER
                DO JL=KIDIA,KFDIA
        
                    ! no on line dry deposition  
                    ZAERO = PAERI(JL,KLEV,JAER) + PTSPHY * PTAERI(JL,KLEV,JAER)
                    ZLAND = PLSM(JL)
                    ZOCEA = 1._JPRB-PLSM(JL)
                    ZFROC = ZOCEA*(1._JPRB-PCI(JL))
                    ZSICE = ZOCEA*PCI(JL)

                    ZVDEP(JAER)= ZFROC * YAEROSOLS(JAER)%RDRYDEPVELOCE + &
                               & ZLAND * YAEROSOLS(JAER)%RDRYDEPVELLND + &
                               & ZSICE * YAEROSOLS(JAER)%RDRYDEPVELSIC
                    ! apply diurnal cycle and stability check 
                    ZVDEP(JAER)= MIN( ZMAXVDRY(JL),  ZVDEP(JAER) )
                    
                    !   SR 03/2018 reduce dry deposition over smoother snow surfaces
                    IF (PSNS(JL) > 1.E-3_JPRB .OR. PCI(JL) > 0.5_JPRB) THEN
                        ZVDEP(JAER) = ZVDEP(JAER)/2.5_JPRB
                        ZVDEP(JAER) = MIN(ZVDEP(JAER),3.E-4_JPRB)
                    ENDIF

                    ! RVDPxxx (from su_aerp, derived from LMDZ in m s-1)
                    !   formula is consistent: (xx m-2 s-1) - (m s-1) * (xx kg-1) * (kg m-3)
                    !   PFAERO is then in xx m-2 s-1
        
                    !PVDEP(JL,JAER) = ZVDEP(JAER)
        
        
                    !-- only the surface flux is diminished of the equivalent effect of the dry deposition
                    !   but the vertical distribution of tendencies is untouched 
        
                    !PFAERO(JL,JAER) = ZVDEP(JAER) * ZAERO * PRHO(JL)
                    YAEROUT_ACLIB%PDRYDEP(JL,JAER+KAER_START-1) = ZVDEP(JAER) * ZAERO * PRHO(JL)
                    PSFAER(JL,JAER) = PSFAER(JL,JAER) + YAEROUT_ACLIB%PDRYDEP(JL,JAER+KAER_START-1)

                    ! print*, JAER, JL, ZVDEP(JAER), ZAERO,  PRHO(JL), ZMAXVDRY(JL), PTSPHY, PSFAER(JL,JAER), PTAERI(JL,KLEV,JAER)

                ENDDO
        
              ENDIF
            ENDIF
        
        ENDDO
        
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_DRYDEP:TACTIC_DRYDEP',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
    END SUBROUTINE TACTIC_DRYDEP

END MODULE M_TACTIC_DRYDEP
