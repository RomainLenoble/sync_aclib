MODULE M_TACTIC_EMIS_SALT

    CONTAINS

    SUBROUTINE TACTIC_EMIS_SALT &
        !---inputs
        & (YDCST_ACLIB, YAEROSOLS, &
        &  KIDIA, KFDIA, KLON, KAERO, KLUOUT, LDSST, &
        &  PCI  , PLSM , PCLK , PWIND, PSST, PRH, &
        !---outputs
        &  PFLXSS)
    
        !*** * FROM AER_SSALT_TACTIC* - SOURCE TERMS FOR SEA SALT AEROSOLS
        
        !**   INTERFACE.
        !     ----------
        !          *TACTIC_EMIS_SALT* IS CALLED FROM *AER_ACLIB_EMIS*.
        
        !     AUTHOR.
        !     -------
        !        Original version (aer_ssalt_tactic)
        !        P. Nabat April 2023
        
        !     SOURCE.
        !     -------
        !     Implementation of the parameterization of Grythe et al. (2014,
        !     page 1286 eq. 7), in which sea salt emission fluxes depend both on
        !     surface wind speed, and possibly SST.
        
        !     MODIFICATIONS.
!-----------------------------------------------------------------------

        USE PARKIND1,   ONLY : JPIM, JPRB
        USE YOMHOOK,    ONLY : LHOOK,DR_HOOK

        USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
        USE MODD_PARAM_AER_ACLIB, ONLY : TAEROSOL

        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        
        !*       0.1  ARGUMENTS
        !             ---------
        
        !---input fields
        ! STRUCTURES
        TYPE(TCST_ACLIB),                INTENT(IN)  :: YDCST_ACLIB
        TYPE(TAEROSOL),DIMENSION(KAERO), INTENT(IN)  :: YAEROSOLS

        INTEGER(KIND=JPIM),              INTENT(IN)  :: KIDIA
        INTEGER(KIND=JPIM),              INTENT(IN)  :: KFDIA
        INTEGER(KIND=JPIM),              INTENT(IN)  :: KLON
        INTEGER(KIND=JPIM),              INTENT(IN)  :: KAERO
        INTEGER(KIND=JPIM),              INTENT(IN)  :: KLUOUT
        
        LOGICAL           ,              INTENT(IN)  :: LDSST  ! !Emissions depending on SST
        REAL(KIND=JPRB)   ,              INTENT(IN)  :: PCI(KLON), PLSM(KLON), PCLK(KLON)
        REAL(KIND=JPRB)   ,              INTENT(IN)  :: PWIND(KLON)
        REAL(KIND=JPRB)   ,              INTENT(IN)  :: PSST(KLON)
        REAL(KIND=JPRB)   ,              INTENT(IN)  :: PRH(KLON)
        
        REAL(KIND=JPRB)   ,              INTENT(OUT) :: PFLXSS(KLON,KAERO)
        
        !-----------------------------------------------------------------------
        
        INTEGER(KIND=JPIM)::ISIZE
        INTEGER(KIND=JPIM) :: JL
        
        REAL(KIND=JPRB)::ZRAYM,ZEMSALT,ZTSC,ZSALTMAS,ZMMD,ZMMDCM,ZRH
        REAL(KIND=JPRB)::ZLIM1,ZLIM2,ZLIM3,ZLIM4       !Bin limits (dry diameter)
        REAL(KIND=JPRB)::ZSS_C1,ZSS_C2,ZSS_C3,ZSS_C4   
        REAL(KIND=JPRB)::ZCOEFSST(KLON)
        REAL(KIND=JPRB)::ZSSC3
        REAL(KIND=JPRB) :: ZFROC, ZOCEA
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_EMIS_SALT:TACTIC_EMIS_SALT',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------

      ! print*, "Tactic emiss salt LSM", MAXVAL(PLSM), MINVAL(PLSM), MAXVAL(PSST), MINVAL(PSST)

        PFLXSS(KIDIA:KFDIA,1:KAERO) = 0._JPRB
        
        ! Dependence on SST
        IF (LDSST) THEN  ! fluxes depend on sst
           DO JL=KIDIA,KFDIA
             ZTSC = PSST(JL)-273.15_JPRB
             ZCOEFSST(JL) = 0.3_JPRB + 0.1_JPRB*ZTSC -0.0076_JPRB*ZTSC**2._JPRB + 0.00021_JPRB*ZTSC**3._JPRB
           ENDDO
        ELSE
          ZCOEFSST(KIDIA:KFDIA)=1.0
        ENDIF
        
        ! Size bin limits
         ZLIM1=0.01_JPRB
         ZLIM2=1.0_JPRB
         ZLIM3=10.0_JPRB
         ZLIM4=100.0_JPRB
        
        ! Parameters for sea-salt emissions
         ZSS_C1=0.7674_JPRB
         ZSS_C2=3.079_JPRB
         ZSS_C3=2.573E-11_JPRB
         ZSS_C4=-1.424_JPRB
        
        ! Sea-salt emissions
        DO JL=KIDIA,KFDIA
          ZSALTMAS=0._JPRB
        
          ! LSM treats lakes as ocean, but (with a few exceptions e.g. the Dead Sea) they don't produce sea salt.
          ! This logic should work for either a fractional LSM or a binary one
          IF (PLSM(JL) == 0._JPRB) THEN
            ! Even with a fractional LSM, this must be all ocean or all lake, because there would be no land to divide them.
            ! Thresholding ensures sensible behaviour if LSM is binary but lake cover fractional.
            IF (PCLK(JL) >= 0.5_JPRB) THEN
              ZOCEA= 0._JPRB
            ELSE
              ZOCEA= 1._JPRB
            ENDIF
          ELSE
            ZOCEA= MAX(1._JPRB-PLSM(JL)-PCLK(JL), 0._JPRB)
          ENDIF
          ! Also exclude sea ice
          ZFROC= ZOCEA*(1._JPRB-PCI(JL))
          !-- flux is considered only over full ocean grids, and for their open ocean
          !  fraction
          IF (ZOCEA == 1._JPRB) THEN
             ZSSC3=ZSS_C3*(1+0.004*(298.15-PSST(JL)))    !c3=2.573E-11
             ZRH = MAX(1E-5,MIN(0.99,PRH(JL)))
        
             ! First bin : film droplets (Grythe et al. 2014)
               ISIZE=1
               ZMMD=YAEROSOLS(ISIZE)%RDIAMED ! in m
               ZMMDCM = ZMMD * 100.  ! in cm
               ZRAYM=(ZSS_C1*(ZMMDCM**ZSS_C2)/(ZSSC3*(ZMMDCM**ZSS_C4)-log10(ZRH)) &
                   & +ZMMDCM**3)**(1./3.)*1.0E4 ! in um (Gerber, 1985)
               ZEMSALT = 235*PWIND(JL)**3.5*exp(-0.55*(log(ZRAYM/0.1))**2) &
                   & + 0.2*PWIND(JL)**3.5*exp(-1.5*(log(ZRAYM/3))**2) &
                   & + 6.8*PWIND(JL)**3*exp(-1.0*((log(ZRAYM/30))**2))
               ZSALTMAS= 4./3.*YDCST_ACLIB%RPI*ZEMSALT*YAEROSOLS(ISIZE)%RDENSITY*ZMMD**3*(ZLIM2-ZLIM1)/2
               PFLXSS(JL,ISIZE) = ZFROC*ZSALTMAS*ZCOEFSST(JL)

             ! Second bin : jet droplets (Grythe et al. 2014)
               ISIZE=2
               ZMMD=YAEROSOLS(ISIZE)%RDIAMED ! in m
               ZMMDCM = ZMMD * 100.  ! in cm
               ZRAYM=(ZSS_C1*(ZMMDCM**ZSS_C2)/(ZSSC3*(ZMMDCM**ZSS_C4)-log10(ZRH)) &
                   & +ZMMDCM**3)**(1./3.)*1.0E4 ! in um
               ZEMSALT = 235*PWIND(JL)**3.5*exp(-0.55*(log(ZRAYM/0.1))**2) &
                   & + 0.2*PWIND(JL)**3.5*exp(-1.5*(log(ZRAYM/3))**2) &
                   & + 6.8*PWIND(JL)**3*exp(-1.0*((log(ZRAYM/30))**2))
               ZSALTMAS= 4./3.*YDCST_ACLIB%RPI*ZEMSALT*YAEROSOLS(ISIZE)%RDENSITY*ZMMD**3*(ZLIM3-ZLIM2)/2
               PFLXSS(JL,ISIZE) = ZFROC*ZSALTMAS*ZCOEFSST(JL)

             ! Third bin : spume droplets (Andreas et al. 1998)
               ISIZE=3
               ZMMD=YAEROSOLS(ISIZE)%RDIAMED ! in m
               ZMMDCM = ZMMD * 100.  ! in cm
               ZRAYM=(ZSS_C1*(ZMMDCM**ZSS_C2)/(ZSSC3*(ZMMDCM**ZSS_C4)-log10(ZRH)) &
                    & +ZMMDCM**3)**(1./3.)*1.0E4 ! in um
               ZEMSALT = 3.5*(10**(0.0676*PWIND(JL)+2.43)*exp(-3.1*(log(ZRAYM/2.1))**2) &
                    & + 10**(0.959*PWIND(JL)**0.5-1.476)*exp(-3.3*(log(ZRAYM/9.2))**2))
               ZSALTMAS= 4./3.*YDCST_ACLIB%RPI*ZEMSALT*YAEROSOLS(ISIZE)%RDENSITY*ZMMD**3*(ZLIM4-ZLIM3)/2
               PFLXSS(JL,ISIZE) = ZFROC*ZSALTMAS*ZCOEFSST(JL)
        
          ENDIF
        ENDDO
        
        
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_EMIS_SALT:TACTIC_EMIS_SALT',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
    END SUBROUTINE TACTIC_EMIS_SALT

END MODULE M_TACTIC_EMIS_SALT
