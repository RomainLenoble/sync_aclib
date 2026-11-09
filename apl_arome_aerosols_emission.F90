SUBROUTINE APL_AROME_AEROSOLS_EMISSION(YDCST, YDMODEL, YDCPG_BNDS, YDCPG_OPTS, & 
                        & YDMF_PHYS_BASE_STATE, &
                        & YDMF_PHYS_SURF, &
                        & YDMF_PHYS, &
                        & YDCPG_MISC, &
                        & YDVARS, &
                        & YDCPG_DYN0, &
                        & KIOFF_MFSHAL, &
                        & PZZ, PDZZ, PPABSM, PRHM, PDT, PQS, &
                        & PTHM, PRHODREFM, PAEROMIN, PSFSV, &
                        & PTENDAERO, &
                        & PAEROM, YLAEROUT_ACLIB)

!**** *APL_AROME_AEROSOLS_EMISSION* - Aerosols emissions

!     Author.
!     -------
!      Rolf Heilemann Myhre *Met Norway*
!      Original : 01-09-2023
!      See APL_AROME for modification history
!
!      03-2025 R. Lenoble : Add ACLIB library and tactic scheme for aerosols
!-----------------------------------------------------------------------

! -   INPUT ARGUMENTS.
!     -------------------
!
!   YDCST:                Model constants
!   YDMODEL:              Model configuration settings
!   YDCPG_BNDS:           Array dimensions and bounds
!   YDCPG_OPTS:           Various dimensions
!   YDMF_PHYS_BASE_STATE: Base state variables (thermodynamic fields)
!   YDMF_PHYS_SURF:       Surface physics variables
!   YDMF_PHYS:            Physics fields
!   YDCPG_DYN0:           Dynamic fields
!   YDCPG_MISC:           Diagnostic variables (for nebulosity) ?
!   YDVARS:               Persistent fields
!   KIOFF_MFSHAL:         Offset for first useful field (EZDIAG use in Chemistry/Dust), (=3)
!   PZZ:                  [m] Height of layers
!   PDZZ:                 [m] Layer thickness
!   PPABSM:               YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF
!   PRHM:                 [kg/m3] Air volumic mass
!   PDT:                  [s] Time step
!   PTHM:                 [K] Potential temperature
!   PRHODREFM:            [kg/m3] Density of air (reference state)
!   PAEROMIN:             [moments/molec_air] Transported moments of dust/aerosols
!
! -   INPUT/OUTPUT ARGUMENTS.
!     -------------------------
!
!   PTENDAERO:            [units/s] Passive scalars tendency (AEROSOL TENDENCY)
!   PSFSV:                [units/m2/s] Passive scalar surface flux for aerosols
!   YLAEROUT_ACLIB        Type ACLIB for aerosols

! -   OUTPUT ARGUMENTS.
!     ------------------
!
!   PAEROM:               [units] Passive scalars (AEROSOL CONCENTRATION)

  USE PARKIND1, ONLY: JPRB, JPIM

  USE YOMCST,                        ONLY: TCST
  USE TYPE_MODEL,                    ONLY: MODEL
  USE CPG_OPTS_TYPE_MOD,             ONLY: CPG_BNDS_TYPE, CPG_OPTS_TYPE
  USE MF_PHYS_BASE_STATE_TYPE_MOD,   ONLY: MF_PHYS_BASE_STATE_TYPE
  USE CPG_TYPE_MOD,                  ONLY: CPG_DYN_TYPE, CPG_MISC_TYPE
  USE MODD_NSV,                      ONLY: NSV_DSTBEG, NSV_DSTEND
  USE MODD_CST_ACLIB,                ONLY: YDCST_ACLIB
  USE MODD_ACLIB_n,                  ONLY: YPAER_ACLIB, YAEROSOLS, YAEROUT_ACLIB
  USE MF_PHYS_SURFACE_TYPE_MOD,      ONLY: MF_PHYS_SURF_TYPE
  USE MODD_ACLIB_N,                  ONLY: YPCHEM_ACLIB, YCHEMS
  USE FIELD_VARIABLES_MOD,           ONLY: FIELD_VARIABLES
  USE MODD_PARAM_CHEM_ACLIB,         ONLY : TPARAM_CHEM_ACLIB
  USE MF_PHYS_TYPE_MOD,              ONLY: MF_PHYS_TYPE
  USE YOMLUN,                        ONLY : NULOUT
  USE MODD_PARAM_AER_ACLIB,          ONLY : TAEROUT_ACLIB
  USE M_INIT_AER_MODEL,              ONLY : INIT_AER_MODEL
  USE MODD_SURFEX_ARO, ONLY : YSURFEX_ARO_CUR, YSURFEX_ARO_ALL
  USE MODI_GET_SURF_VAR_N
  USE MODD_FROMMPA,      ONLY : CPROGRAM_MSE

  USE YOMHOOK,                  ONLY: LHOOK, DR_HOOK, JPHOOK

  IMPLICIT NONE

  TYPE(TCST),                    INTENT(IN)    :: YDCST
  TYPE(MODEL),                   INTENT(IN)    :: YDMODEL
  TYPE(CPG_BNDS_TYPE),           INTENT(IN)    :: YDCPG_BNDS
  TYPE(CPG_OPTS_TYPE),           INTENT(IN)    :: YDCPG_OPTS
  TYPE(MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
  TYPE(MF_PHYS_SURF_TYPE),       INTENT(IN)    :: YDMF_PHYS_SURF
  TYPE(CPG_DYN_TYPE),            INTENT(IN)    :: YDCPG_DYN0
  TYPE(CPG_MISC_TYPE),           INTENT(IN)    :: YDCPG_MISC
  TYPE(FIELD_VARIABLES),         INTENT(IN)    :: YDVARS
  TYPE(MF_PHYS_TYPE),            INTENT(IN)    :: YDMF_PHYS

  INTEGER(KIND=JPIM), INTENT(IN) :: KIOFF_MFSHAL

  REAL(KIND=JPRB), INTENT(IN)    :: PZZ(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB), INTENT(IN)    :: PDZZ(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB), INTENT(IN)    :: PPABSM(YDCPG_BNDS%KFDIA, 0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB), INTENT(IN)    :: PRHM(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB), INTENT(IN)    :: PDT
  REAL(KIND=JPRB), INTENT(IN)    :: PQS(YDCPG_BNDS%KFDIA)       ! surface humidity

  REAL(KIND=JPRB), INTENT(IN)    :: PTHM(YDCPG_BNDS%KFDIA, 0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB), INTENT(IN)    :: PRHODREFM(YDCPG_BNDS%KFDIA, 0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB), INTENT(IN)    :: PAEROMIN(YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NAERO)

  REAL(KIND=JPRB), INTENT(INOUT) :: PSFSV(YDCPG_BNDS%KFDIA, YDMODEL%YRML_GCONF%YGFL%NAERO)

  REAL(KIND=JPRB), INTENT(INOUT) :: PTENDAERO(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NAERO)
  TYPE(TAEROUT_ACLIB), INTENT(INOUT) :: YLAEROUT_ACLIB

  REAL(KIND=JPRB), INTENT(OUT)   :: PAEROM(YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NAERO)


  INTEGER(KIND=JPIM) :: IEZDIAG_CHEM
  REAL(KIND=JPRB) :: ZPEZDIAG(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG - KIOFF_MFSHAL + 1)

  INTEGER(KIND=JPIM) :: JLON, JLEV, JGFL, KBL, JI
  INTEGER(KIND=JPIM), PARAMETER :: IKL = -1

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-- aerosol-related quantities for tactic
  INTEGER(KIND=JPIM),PARAMETER    :: JPTILES=8_JPIM
  REAL(KIND=JPRB)   :: ZCI(YDCPG_BNDS%KFDIA) ! sea ice fraction in ecmwf physics, but here sea ice
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)                               :: ZAERLTS, ZAERSCC 
  REAL(KIND=JPRB)   :: ZWND(YDCPG_BNDS%KFDIA)
  REAL(KIND=JPRB),DIMENSION(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)           :: ZTP   ! Temperature for ACLIB

  ! Surfex variables
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: PGZ0
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: PBARESOIL
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: ZCLK 
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: ZSEA
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: ZNATURE
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: PLSM
  REAL(KIND=JPRB), DIMENSION(YDCPG_BNDS%KFDIA)    :: PTDSNOW   ! Snow depth

  INTEGER, ALLOCATABLE, DIMENSION(:) :: KMASK_SFX

  !------------------------------------------------------------------

#include "aro_mnhdust.h"
#include "aer_main1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('APL_AROME_AEROSOLS_EMISSION_MOD:APL_AROME_AEROSOLS_EMISSION', 0, ZHOOK_HANDLE)

  ASSOCIATE(YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1)
  ASSOCIATE(TMERGL=>YDPHY1%TMERGL, GCONV=>YDPHY1%GCONV, RD1=>YDPHY1%RD1, NDDUST=>YDMODEL%YRML_PHY_AER%YREAERSRC%NDDUST,RG=>YDCST%RG,NAERO=>YDMODEL%YRML_GCONF%YGFL%NAERO)
  
  IEZDIAG_CHEM = YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG-KIOFF_MFSHAL + 1
  PAEROM(:,:,:) = PAEROMIN(:,:,:)


  ! Configuration with tactic
  IF (TRIM(YDMODEL%YRML_CHEM%YRCOMPO%AERO_SCHEME)=="tactic") THEN
    
    KBL = YDCPG_BNDS%KBL
    YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (KBL)

    ZCI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = 0.0_JPRB   ! Ice fraction
    ZCLK(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = 0.0_JPRB  ! Lake fraction

    CALL GET_SURF_VAR_N &
    (YSURFEX_ARO_CUR%FM, YSURFEX_ARO_CUR%IM,    &
    YSURFEX_ARO_CUR%SM, YSURFEX_ARO_CUR%TM, &
    YSURFEX_ARO_CUR%WM, YSURFEX_ARO_CUR%DUO,     &
    YSURFEX_ARO_CUR%DU, YSURFEX_ARO_CUR%UG,  &
    YSURFEX_ARO_CUR%U,      YSURFEX_ARO_CUR%USS,     &
    YSURFEX_ARO_CUR%SM%S,                            &
    CPROGRAM_MSE, YDCPG_BNDS%KFDIA,1, &
        PZ0=PGZ0,                               &  ! Roughness length
        PBARE=PBARESOIL,                        &  ! Bare soil fraction
        PSEA=ZSEA,                              &  ! Land sea mask
        PNATURE=ZNATURE,                        &  ! Patch nature (necessaire pour la routine)
        PTDSNOW=PTDSNOW)                           ! Snow depth


    ! get sand and clay fraction
    KMASK_SFX = YSURFEX_ARO_CUR%U%NR_NATURE
    DO JI=1, YSURFEX_ARO_CUR%U%NSIZE_NATURE
      JLON=KMASK_SFX(JI) 
      ZAERSCC(JLON)=YSURFEX_ARO_CUR%IM%K%XCLAY(JI, 1)*100
      ZAERLTS(JLON)=YSURFEX_ARO_CUR%IM%K%XSAND(JI, 1)*100
    ENDDO

    ! convert LSM from 1=sea to 0=sea
    PLSM = 1.-ZSEA
    ! PLSM = ZSEA

    ! Wind speed at 10m
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZWND(JLON) = SQRT(YDMF_PHYS%OUT%UCLS(JLON)**2 + YDMF_PHYS%OUT%VCLS(JLON)**2)
    ENDDO

  ! Temperature computation for AER_ACLIB
    ZTP(:,:)=0.
    ZTP(:,:)=PTHM(:,:)*(PPABSM(:,:) / YDCST_ACLIB%RP00)**(YDCST_ACLIB%RRD/YDCST_ACLIB%RCPD)


    CALL INIT_AER_MODEL(YPAER_ACLIB, &
    & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NULOUT, &
    & YLAEROUT_ACLIB)

    ! Tactic emissions are negative compare to arome surface. Need to change sign before/after tactic
    PSFSV(:,:) = PSFSV(:,:) * -1.

    CALL AER_MAIN1( YDCST_ACLIB=YDCST_ACLIB, YPAER_ACLIB=YPAER_ACLIB, &
      & YAEROSOLS=YAEROSOLS, &
      & KIDIA=YDCPG_BNDS%KIDIA, & 
      & KFDIA=YDCPG_BNDS%KFDIA, &
      & KLUOUT=NULOUT, &
      & KLON=YDCPG_OPTS%KLON, KLEV=YDCPG_OPTS%KFLEVG, PTSTEP=YDCPG_OPTS%ZDTPHY, &
      & PAP=PPABSM(:, 1:YDCPG_OPTS%KFLEVG), &                                              ! FULL-LEVEL PRESSURE
      & PAPH=YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:), &     ! HALF-LEVEL PRESSURE    
      & PAPHI=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:), &    ! GEOPOTENTIAL ON HALF-LEVELS
      & PDELP=YDCPG_DYN0%XYB%DELP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:), &                  ! DELTA PRESSURE
      & PTP=ZTP,  &                 ! FULL-LEVEL TEMPERATURE
      & PVERVEL=YDMF_PHYS_BASE_STATE%YCPG_PHY%W(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:),&     ! FULL-LEVEL VERTICAL VELOCITY
      & PRH=PRHM, &                                   ! RELATIVE HUMIDITY
      & PRHODREFM=PRHODREFM, &                        ! Dry air density
      & PSST=YDMF_PHYS_BASE_STATE%YGSP_RR%T, &        ! Sea surface temperature
      & PWS1=PQS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),& ! Soil humidity, need to compute soil humidity
      & PCI=ZCI, &                                    ! FRACTION OF SEA-ICE, need to compute fraction of sea ice
      & PBARESOIL= PBARESOIL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),&                 ! BARE SOIL FRACTION
      & PLSM=PLSM,&      ! LAND-SEA MASK
      & PCLK=ZCLK(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &             ! Water fraction, equal to 0 ?
      & PSNS=PTDSNOW(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &          ! SNOW DEPTH
      & PAERDEP=ZAERLTS, &                                          ! SAND CONTENT
      & PAERSCC=ZAERSCC, &                                          ! CLAY CONTENT
      & PZ0M=PGZ0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &             ! ROUGHNESS LENGTH
      & PWIND=ZWND,&                                                ! SURFACE WIND SPEED           (m/s)
      & PAER=PAEROM,&                                               ! AEROSOL CONCENTRATION
      & PTENAER=PTENDAERO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,:),&  ! AEROSOL TENDENCY
      & PSFAER=PSFSV,&                                              ! AEROSOL SURFACE FLUX
      & YAEROUT_ACLIB=YLAEROUT_ACLIB)                               ! STRUCTURE FOR AEROSOL OUTPUT

  ENDIF

  ! RL: need to have positive values
  PSFSV(:,:) = PSFSV(:,:) * -1.

  END ASSOCIATE
  END ASSOCIATE

  IF (LHOOK) CALL DR_HOOK('APL_AROME_AEROSOLS_EMISSION_MOD:APL_AROME_AEROSOLS_EMISSION', 1, ZHOOK_HANDLE)

END SUBROUTINE APL_AROME_AEROSOLS_EMISSION
