SUBROUTINE ARO_GROUND_PARAM (KBL, KLON, KST, KEND, KSTEP, KRR, KSW, KSV, &
                           & KDGUNG, KDGUXG, KDLUNG, KDLUXG, &
                           & LD_MPA, CCOUPLING,  &
                           & LDXFUMSE, KINDAT, PRHGMT, PTIMEC, PSOVR, PCOSDEL, PSINDEL,   &
                           & PINDX, PINDY, PUM, PVM, PTM, PQVM, PSVM, PCCO2, PRHODREF,    &
                           & PPABSM, PSURF, PTSTEP, PZZ, PZS, PZSEPS, PMU0,PMU0N, PLON,   &
                           & PGEMU, PSW_BANDS, PINPRR, PINPRS, PINPRG, PFRTHDS, PFRSOSCDS,&
                           & PFRSODIDS, PCFAQ, PCFATH, PCFAU, PCFBQ,PCFBTH, PCFBU, PCFBV, &
                           & PSFTH, PSFRV, PSFSV, PSFCO2, PSFU, PSFV, PDIR_ALB, PSCA_ALB, &
                           & PEMIS, PTSRAD, LWNETF,YDMODEL)

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!****  *AC_GROUND_PARAM*

!    PURPOSE
!    -------
!       Monitor to call the externalized surface from AROME (adaptation of
!       GROUND_PARAM of Meso-NH)

!**  METHOD
!    ------

!    EXTERNAL
!    --------

!    IMPLICIT ARGUMENTS
!    ------------------

!    REFERENCE
!    ---------
!      Noilhan and Planton (1989)

!    AUTHOR
!    ------
!       S. Malardel           * Meteo-France *

!    MODIFICATIONS
!    -------------
!       Original    07/03 (1D AROME prototype)
!       phasing for 3D AROME prototype 03/04
!       R. Zaaboul 15-Feb-2006 add new arguments KBL,KGPCOMP and LD_NBLOCKTOT and modification of NPROMA loop
!       L. Kraljevic 09/06 add new arguments PCFAQ, PCFATH, PCFAU, PCFBQ,PCFBTH, PCFBU, PCFBV and modifications for ALADIN
!       L. Kraljevic 09/06 removed reshapings of arrays
!       A. ALias 01-2007   :
!            - Add test on KSV
!            - Do not go through LUSECHEM_MSE/LCH_SURFACE_FLUX_MSE as long as
!              they are not initialised
!       Y. Seity   05/07 remove chemical fields conversions
!       F. Bouyssel 05/07 Implicit coupling for wind when CCOUPLING=E + Cleaning
!       F. Bouyssel 03/08 Cleaning for CCOUPLING=E,I,V
!       Y. Seity   09/08 modification for DUSTS
!       R. El Khatib 06/10 Resources savings
!       P.Marguinaud 08/10 OpenMP + cleaning
!       F. Bouyssel 05/10 Bug correction and use of PINDX, PINDY
!       B. Decharme 12/10 bugfix Qair <0.0 and PMU0N added
!       F. Bouyssel 05/11 Initialization of ZCFAQ
!       A. Alias    05/11 threshold added for computation of ZQA
!       P. Marguinaud 09/12 control threshold for orography; check SURFEX orography vs model
!       C. Wastl 09/14 Adaptations for orographic shadowing
!       J. Masek 08/2023 Reintroduced orography check, but with removed useless
!                        call to MPL_BARRIER causing deadlock (not activated)
!       M. Cussac, M. Michou : Modifications to consider aerosol and/or chemical scheme of ARPEGE-Climat

!-------------------------------------------------------------------------------

!*       0.     DECLARATIONS
!               ------------

USE MODD_SURFEX_ARO, ONLY : YSURFEX_ARO_CUR, YSURFEX_ARO_ALL
USE MODD_CSTS,          ONLY : XP00, XCPD, XRD, XRV, XRHOLW, XDAY, XPI, XLVTT, XMD, XAVOGADRO, XSTEFAN

USE MODD_SURF_PAR,     ONLY :  XUNDEF

USE MODD_IO_SURF_ARO,  ONLY : XLON0,XBETA,XRPK,NLUOUT
USE MODD_FROMMPA,      ONLY : LCARTESIAN_MSE, CPROGRAM_MSE,         &
                              CSV_MSE, NSV_CHEMBEG_MSE, NSV_CHEMEND_MSE,&
                              NSV_AERBEG_MSE, NSV_AEREND_MSE,       &
                              NSV_DSTBEG_MSE, NSV_DSTEND_MSE,       &
                              NSV_DSTDEPBEG_MSE, NSV_DSTDEPEND_MSE, &
                              NSV_CO2_MSE,                      &
                              LDUST_MSE,LDEPOS_DST_MSE,LORILAM_MSE,     &
                              LUSECHEM_MSE,              &
                              LCH_SURFACE_FLUX_MSE
USE MODD_DST_SURF,     ONLY : XMOLARWEIGHT_DST
USE MODI_GET_SURF_VAR_N
USE MODI_CORRECT_TIME_FLAKE
USE MODI_CORRECT_TIME_ISBA
USE MODI_CORRECT_TIME_SEA
USE MODI_CORRECT_TIME_SURF
USE MODI_CORRECT_TIME_TEB
USE MODI_CORRECT_TIME_WAT
USE MODI_COUPLING_SURF_ATM_N

USE MODD_ISBA_N
USE MODD_PREP_SNOW, ONLY : NIMPUR


USE YOMCT0, ONLY : LELAM
USE TYPE_MODEL     , ONLY : MODEL

IMPLICIT NONE

!*      0.1    declarations of arguments

INTEGER(KIND=JPIM),          INTENT(IN):: KBL       ! NPROMA-packet number
INTEGER(KIND=JPIM),          INTENT(IN):: KLON      ! horizontal dimension
INTEGER(KIND=JPIM),          INTENT(IN):: KST       ! start of physical values
INTEGER(KIND=JPIM),          INTENT(IN):: KEND      ! end of physical values
INTEGER(KIND=JPIM),          INTENT(IN):: KSTEP
INTEGER(KIND=JPIM),          INTENT(IN):: KRR       ! number of hydrometeors
INTEGER(KIND=JPIM),          INTENT(IN):: KSW       ! number of solar spectral band
INTEGER(KIND=JPIM),          INTENT(IN):: KSV       ! number of passive scalar
INTEGER(KIND=JPIM),          INTENT(IN):: KDGUNG, KDGUXG, KDLUNG, KDLUXG ! dim param
                                         ! for E-zone
LOGICAL,          INTENT(IN):: LD_MPA    ! LMPA from yomarphy
CHARACTER (LEN=1),INTENT(IN):: CCOUPLING ! coupling mode 'E' or 'I' or 'V'
LOGICAL,          INTENT(IN):: LDXFUMSE  ! T if CDCONF(1)=X (XFU 0 time step)

!* date and angles

INTEGER(KIND=JPIM),          INTENT(IN):: KINDAT    ! date and time of model start in ALADIN form
REAL(KIND=JPRB),             INTENT(IN):: PRHGMT    ! time in seconds from 00TU
REAL(KIND=JPRB),             INTENT(IN):: PTIMEC    !  current time since midnight (UTC, s)
REAL(KIND=JPRB),             INTENT(IN):: PSOVR     ! true solar time, at longiture 0, in seconds
REAL(KIND=JPRB),             INTENT(IN):: PCOSDEL   ! cosine of the solar declination angle
REAL(KIND=JPRB),             INTENT(IN):: PSINDEL   ! sine of the solar declination angle

!* Input fields
!  ------------
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PINDX     ! x coordinate on the plane
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PINDY     ! y coordinate on the plane
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PUM,PVM   ! wind components at the instant t-dt
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PTM       ! theta
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PQVM      ! water mixing ratio
REAL(KIND=JPRB), DIMENSION(KLON,KSV),INTENT(IN)   :: PSVM      ! passive scalar
REAL(KIND=JPRB),                     INTENT(IN)   :: PCCO2     ! CO2 concentration
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PRHODREF  ! dry air density
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PPABSM    ! pressure
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PSURF     ! surface pressure
REAL(KIND=JPRB),                     INTENT(IN)   :: PTSTEP    ! time-step of atmospheric model
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PZZ       ! height of model intercouches
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PZS       ! orography (of atmospheric model)
REAL(KIND=JPRB),                     INTENT(IN)   :: PZSEPS    ! relative difference between model and SURFEX orography
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PMU0      ! cos of zenithal angle at t
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PMU0N     ! cos of zenithal angle at t+dt
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PLON      ! longitude in radian (!! in MesoNH, angle are usually in degre)
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PGEMU     ! sin of latitude
REAL(KIND=JPRB), DIMENSION(KSW),     INTENT(IN)   :: PSW_BANDS ! solar spectral bands
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PINPRR    ! liquid precipitation surface flux
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PINPRS    ! snow precipitation surface flux
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PINPRG    ! graupel precipitation surface flux
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(INOUT)   :: PFRTHDS   ! flux thermique descendant a la surface
REAL(KIND=JPRB), DIMENSION(KLON,KSW),INTENT(INOUT)   :: PFRSOSCDS ! flux solaire diffus descendant a la surface for each spectral band
REAL(KIND=JPRB), DIMENSION(KLON,KSW),INTENT(INOUT)   :: PFRSODIDS ! flux solaire direct descendant a la surface for each spectral band
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFAQ
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFATH
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFAU
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFBQ
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFBTH
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFBU
REAL(KIND=JPRB), DIMENSION(KLON),    INTENT(IN)   :: PCFBV
TYPE(MODEL)    ,                     INTENT(IN)   :: YDMODEL

!* Output fields
!  -------------

REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(OUT) :: PSFTH     ! surface flux of potential temperature (Km/s)
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(OUT) :: PSFRV     ! surface flux of water vapor           (m/s*kg/kg)
REAL(KIND=JPRB), DIMENSION(KLON,KSV), INTENT(OUT) :: PSFSV     ! surface flux of scalar                (m/s*kg/kg)
                                                      ! flux of chemical var.                 (ppp.m/s)
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(OUT) :: PSFCO2    ! surface flux of CO2                   (m/s*kg/kg)
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(OUT) :: PSFU      ! surface fluxes of horizontal
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(OUT) :: PSFV      ! momentum in x and y directions        (m2/s2)

REAL(KIND=JPRB), DIMENSION(KLON,KSW), INTENT(INOUT) :: PDIR_ALB  ! direct  albedo for each spectral band (-)
REAL(KIND=JPRB), DIMENSION(KLON,KSW), INTENT(INOUT) :: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(INOUT) :: PEMIS     ! surface emissivity                    (-)
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(INOUT) :: PTSRAD    ! surface radiative temperature         (K)
!orographic shadowing
REAL(KIND=JPRB), DIMENSION(KLON),     INTENT(INOUT) :: LWNETF    ! longwave net flux            (W/m2)
!-------------------------------------------------------------------------------

!*      0.2    declarations of local variables
!              -------------------------------

!* Atmospheric variables
!  ---------------------

!            suffix 'A' stands for atmospheric variable at first model level

REAL, DIMENSION(KLON)    :: ZZREF   ! Forcing height
REAL, DIMENSION(KLON)    :: ZTA     ! Temperature
REAL, DIMENSION(KLON)    :: ZQA     ! humidity (kg/m3)
REAL, DIMENSION(KLON)    :: ZPA     ! Pressure
REAL, DIMENSION(KLON)    :: ZPS     ! Pressure
REAL, DIMENSION(KLON)    :: ZRAIN   ! liquid precipitation (kg/m2/s)
REAL, DIMENSION(KLON)    :: ZSNOW   ! solid precipitation  (kg/m2/s)
REAL, DIMENSION(KLON)    :: ZTSUN   ! solar time           (s since midnight)

REAL, DIMENSION(KLON)    :: ZUA     ! u component of the wind parallel to the orography
REAL, DIMENSION(KLON)    :: ZVA     ! v component of the wind parallel to the orography
REAL, DIMENSION(KLON)    :: ZU      ! zonal wind
REAL, DIMENSION(KLON)    :: ZV      ! meridian wind
REAL, DIMENSION(KLON)    :: ZWIND   ! wind parallel to the orography
REAL, DIMENSION(KLON)    :: ZRHOA   ! air density
REAL, DIMENSION(KLON)    :: ZDIR    ! wind direction (rad from N clockwise)
REAL, DIMENSION(KLON)    :: ZSFU    ! zonal momentum flux
REAL, DIMENSION(KLON)    :: ZSFV    ! meridian momentum flux
REAL, DIMENSION(KLON)    :: ZCO2    ! CO2 concentration (kg/kg)
REAL, DIMENSION(KLON)    :: ZALFA   ! angle between the wind and the x axis

REAL, DIMENSION(KLON)    :: ZSOLANG ! solar angle
REAL, DIMENSION(KLON)    :: ZZENITH ! zenith angle at t(in rad)
REAL, DIMENSION(KLON)    :: ZZENITH2! zenith angle at t+1 (in rad)
REAL, DIMENSION(KLON)    :: ZSINZEN ! sine of zenith
REAL, DIMENSION(KLON)    :: ZSINAZI ! sine of azimuth
REAL, DIMENSION(KLON)    :: ZCOSAZI ! cosine of azimuth
REAL, DIMENSION(KLON)    :: ZAZIM   ! azimuth (in rad)

REAL, DIMENSION(KLON)    :: ZSFTH   ! Turbulent flux of heat
REAL, DIMENSION(KLON)    :: ZSFTQ   ! Turbulent flux of water
REAL, DIMENSION(KLON)    :: ZSFCO2  ! Turbulent flux of CO2
REAL, DIMENSION(KLON,KSV):: ZSFTS   ! Turbulent flux of scalar
REAL, DIMENSION(KSV)     :: ZMOLMASS  ! passive scalar Molar mass (g/mol)
REAL, DIMENSION(KLON)    :: ZCFAQ
REAL, DIMENSION(KLON)    :: ZCFBQ
REAL, DIMENSION(KLON)    :: ZCFBUV
REAL, DIMENSION(KLON)    :: ZZ1
REAL, DIMENSION(KLON)    :: ZZ2
REAL, DIMENSION(KLON)    :: ZSFMOD
REAL, DIMENSION(KLON,KSV):: ZSVM
REAL, DIMENSION(KLON)    :: ZTSURF    ! surface effective tempe
REAL, DIMENSION(KLON)    :: ZZ0       ! roughness length for mo
REAL, DIMENSION(KLON)    :: ZZ0H      ! roughness length for he
REAL, DIMENSION(KLON)    :: ZQSURF    ! specific humidity at su


CHARACTER(LEN=16),DIMENSION(:), ALLOCATABLE      :: YSV     ! used for chemistry

!* Dimensions
!  ----------

INTEGER :: JLAYER        ! loop counter
INTEGER :: JSV,JI,II     ! loop counter
INTEGER :: ILEN          ! horizontal dimension in usage for Surfex (ie : no useless extra-points)

!* Other integer variables
!  -----------------------

INTEGER :: IYEAR         ! year at initial run time
INTEGER :: IMONTH        ! month at initial run time
INTEGER :: IDAY          ! day at initial run time

!* Other real variables
!  -----------------------
REAL :: ZDGUN            ! y dimension of C+I+E zone
REAL :: ZDGUX            ! y dimension of C+I zone
REAL :: ZDLUN            ! x dimension of C+I+E zone
REAL :: ZDLUX            ! x dimension of C+I zone

REAL, DIMENSION(KLON)    :: ZEMASK  ! E-zone mask
REAL, DIMENSION(KLON,NIMPUR)  :: ZIMPWET      ! 
REAL, DIMENSION(KLON,NIMPUR)  :: ZIMPDRY    

LOGICAL :: LLCHEM_ARPC, LLAERO_ARPC, LLCHEMAERO_ARPC
!-------------------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('ARO_GROUND_PARAM',0,ZHOOK_HANDLE)

ASSOCIATE(LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM, CHEM_SCHEME=>YDMODEL%YRML_CHEM%YRCHEM%CHEM_SCHEME, &
& AERO_SCHEME=>YDMODEL%YRML_CHEM%YRCOMPO%AERO_SCHEME,NAERO=>YDMODEL%YRML_GCONF%YGFL%NAERO)

YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (KBL)

!IF ((KSTEP == 0) .AND. (PZSEPS >= 0)) THEN
!
!  IF (ANY ((ABS (PZS (KST:KEND) - YSURFEX_ARO_CUR%U%XZS (KST:KEND)) > &
!         &  MAX(PZSEPS * ABS (PZS (KST:KEND)),10._JPRB*PZSEPS))             &
!         & .AND. (YSURFEX_ARO_CUR%U%XZS (KST:KEND) > 0._JPRB))) THEN
!
!    WRITE(NLUOUT, '(''ARO_GROUND_PARAM: OROGRAPHY MISMATCH'')')
!    WRITE(NLUOUT, '(5X,2A12)') 'MODEL', 'SURFEX'
!    II = 0
!    DO JI = KST, KEND
!      IF ((ABS (PZS (JI) - YSURFEX_ARO_CUR%U%XZS (JI)) > PZSEPS * ABS (PZS (JI))) &
!         & .AND. (YSURFEX_ARO_CUR%U%XZS (JI) > 0._JPRB)) THEN
!        WRITE(NLUOUT, '(I4,X,2E12.5)') JI, PZS (JI), YSURFEX_ARO_CUR%U%XZS (JI)
!        II = II + 1
!        IF (II > 10) EXIT
!      ENDIF
!    ENDDO
!    CALL ABOR1 ('ARO_GROUND_PARAM: OROGRAPHY MISMATCH')
!  ENDIF
!
!ENDIF


PSFTH    = 0.
PSFRV    = 0.
PSFSV    = 0.
PSFCO2   = 0.
PSFU     = 0.
PSFV     = 0.
!PDIR_ALB = 0.
!PSCA_ALB = 0.
PEMIS    = 0.
PTSRAD   = 0.
ZSVM     = 0.

ZUA      = 0.
ZVA      = 0.
ZALFA    = 0.
ZDIR     = 0.
ZWIND    = 0.
ZU       = 0.
ZV       = 0.
ZPA      = 0.
ZPS      = 0.
ZQA      = 0.
ZRHOA    = 0.
ZRAIN    = 0.
ZSNOW    = 0.
ZTSUN    = 0.
ZZENITH  = 0.
ZZENITH2 = 0.
ZSINZEN  = 0.
ZZREF    = 0.
ZCO2     = 0.

ALLOCATE (YSV(KSV))

! cases of ARPEGE-Climat tropospheric chemical schemes or aerosol scheme 
LLAERO_ARPC     = .FALSE.
LLCHEM_ARPC     = .FALSE.
LLCHEMAERO_ARPC = .FALSE.
LLAERO_ARPC     = AERO_SCHEME == "tactic" .AND. NAERO>0
LLCHEM_ARPC     = LCHEM_ARPCLIM.AND.(TRIM(CHEM_SCHEME).EQ.'arpclim_relacs'.OR.TRIM(CHEM_SCHEME).EQ.'mocage')
LLCHEMAERO_ARPC = LLAERO_ARPC.AND.LLCHEM_ARPC

IF (KSV > 0) THEN
  YSV(:)(1:16)=CSV_MSE(:)(1:16)
ENDIF
!*       1.     CONVERSION OF THE ATMOSPHERIC VARIABLES
!               ---------------------------------------

!        1.2    Horizontal wind direction (rad from N clockwise)
!               -------------------------

ZUA(KST:KEND) = PUM(KST:KEND)
ZVA(KST:KEND) = PVM(KST:KEND)
ZWIND(KST:KEND) = SQRT( ZUA(KST:KEND)**2 + ZVA(KST:KEND)**2 )

!* angle between Y axis and wind (rad., clockwise)

WHERE(ZUA(:)/=0. .OR. ZVA(:)/=0.)
  ZALFA(:)=ATAN2(ZUA(:),ZVA(:))
END WHERE
WHERE(ZALFA(:)<0.) ZALFA(:) = ZALFA(:) + 2. * XPI

!* angle between North and wind (rad., clockwise)

IF (.NOT. LCARTESIAN_MSE) THEN
  ZZ1 (KST:KEND) = PLON(KST:KEND)*(180./XPI)-XLON0
  ZZ2 (KST:KEND) = ZZ1(KST:KEND)-360.*MAX(0.,SIGN(1.,ZZ1(KST:KEND)-180.))
  ZDIR(KST:KEND) = (XRPK*ZZ2(KST:KEND) - XBETA ) * XPI/180. + ZALFA(KST:KEND)
ELSE
  ZDIR(KST:KEND) = - XBETA * XPI/180. + ZALFA(KST:KEND)
ENDIF

ZU(KST:KEND) = ZWIND(KST:KEND) * SIN(ZDIR(KST:KEND))
ZV(KST:KEND) = ZWIND(KST:KEND) * COS(ZDIR(KST:KEND))
!ZU(KST:KEND) = PUM(KST:KEND)
!ZV(KST:KEND) = PVM(KST:KEND)

!        1.8    Temperature from the potential temperature
!               ------------------------------------------

ZTA(KST:KEND) = PTM(KST:KEND)

!        1.6    Pressure and Exner function
!               ---------------------------

ZPS(KST:KEND) = PSURF (KST:KEND)
ZPA(KST:KEND) = PPABSM(KST:KEND)

!        1.7    humidity in kg/m3 from the mixing ratio
!               ---------------------------------------

IF (KRR>0) THEN
  ZQA(KST:KEND) = MAX(0.0,PQVM(KST:KEND) * PRHODREF(KST:KEND))
ENDIF

!        1.9    Air density
!               -----------

ZRHOA(KST:KEND) = PRHODREF(KST:KEND)

!        1.10   Precipitations
!               --------------

IF (.NOT.LD_MPA) THEN
  ZRAIN(KST:KEND) = PINPRR(KST:KEND)
  ZSNOW(KST:KEND) = PINPRS(KST:KEND)
ELSE
  IF (KRR > 2) THEN
    ZRAIN(KST:KEND) = ZRAIN(KST:KEND) + PINPRR(KST:KEND) * XRHOLW
    !!!! Commentaires provisoires, en attendant les precip de la convection
    !IF (CDCONV == 'KAFR')  THEN
    !  ZRAIN = ZRAIN + (PPRCONV - XPRSCONV)  * XRHOLW
    !  ZSNOW = ZSNOW + PPRSCONV * XRHOLW
    !ENDIF
    IF (KRR >= 5) THEN
      ZSNOW(KST:KEND) = ZSNOW(KST:KEND) + PINPRS(KST:KEND) * XRHOLW
    ENDIF
    IF (KRR >= 6) THEN
      ZSNOW(KST:KEND) = ZSNOW(KST:KEND) + PINPRG(KST:KEND) * XRHOLW
    ENDIF
    !!* Schema de grele pas encore implante dans AROME
    !IF(KRR >= 7) THEN
    !  ZSNOW = ZSNOW + XINPRH * XRHOLW
    !ENDIF
  ENDIF
ENDIF

!        1.11.a   Solar time
!               ----------

IF (.NOT. LCARTESIAN_MSE) THEN
   ZTSUN(KST:KEND) = MOD(PSOVR +PLON(KST:KEND)* (180./XPI) *240., XDAY)
ELSE
   ZTSUN(KST:KEND) = MOD(PSOVR +XLON0            *240., XDAY)
ENDIF

ZSOLANG(KST:KEND)=  ((ZTSUN(KST:KEND)/3600.)-12.0)*15.0*(XPI/180.)          ! hour angle in radians

!        1.11.b Zenithal angle
!               --------------

ZZENITH (KST:KEND) = ACOS( PMU0 (KST:KEND) )
ZZENITH2(KST:KEND) = ACOS( PMU0N(KST:KEND) )
ZSINZEN (KST:KEND) = SQRT(1.-PMU0(KST:KEND)**2)


!        1.11.c Azimuthal angle
!               --------------

WHERE (ZSINZEN(:)/=0.)
  ZSINAZI(:) = - PCOSDEL * SIN(ZSOLANG(:)) / ZSINZEN(:)
  ZCOSAZI(:) = (-PGEMU(:)*PCOSDEL*COS(ZSOLANG(:))     &
                   +SQRT(1.-PGEMU(:)**2)*PSINDEL ) / ZSINZEN(:)
  ZAZIM(:) = XPI - ATAN2(ZSINAZI(:),ZCOSAZI(:))
ELSEWHERE
  ZAZIM(:) = 0.
ENDWHERE

!        1.11.d Date
!               ----

!Attention, pour le moment, on passe a la surface uniquement la date initiale
!du run (la seule disponible actuellement dans Aladin).
!Si on veut passer la vraie date a chaque pas de temps, il faut installer
!une routine capable de faire avancer la date. Il en existe une dans mesonh.
!
! ANNEE / MOIS / JOUR DU MOIS
!
IYEAR = KINDAT / 10000
IMONTH = (KINDAT - 10000*IYEAR ) / 100
IDAY = KINDAT - 10000*IYEAR - 100*IMONTH

!        1.12   Forcing level
!               -------------

ZZREF(KST:KEND) = PZZ(KST:KEND)

!        1.13   CO2 concentration (kg/m3)
!               -----------------

ZCO2=0.
IF (NSV_CO2_MSE /= 0) THEN
  ZCO2(KST:KEND) = 1E-6 * PSVM(KST:KEND,NSV_CO2_MSE) * PRHODREF(KST:KEND) * 44./28.8
ELSE
  ZCO2(KST:KEND) = PCCO2 * PRHODREF(KST:KEND)
END IF


!chemical conversion : from part/part to molec./m3

! ZSVM used later in the SURFEX code to compute dry deposition velocity fluxes
IF (KSV > 0) THEN
  !test on chemical scheme in case of arpege climate simulation
  !intructions below should be changed (1) in case of use of dry deposition fluxes of aerosols computed from surfex
  !and (2) in case of dry deposition of both chemistry and aerosols are considered in the same simulation as NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
  !add both aerosols and gazeous species   
  IF (LCHEM_ARPCLIM.AND.TRIM(CHEM_SCHEME).EQ.'arpclim_relacs') THEN
    DO JLAYER=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
      !input is mol/mol
      ZSVM(KST:KEND,JLAYER) = PSVM(KST:KEND,JLAYER) * XAVOGADRO * PRHODREF(KST:KEND) / XMD
    END DO
  ELSEIF (LCHEM_ARPCLIM.AND.TRIM(CHEM_SCHEME).EQ.'mocage') THEN
    DO JLAYER=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
      !input is kg/kg (factor 1000 to account to ZMOLMASS in g/mol)
      ZMOLMASS(JLAYER)=YDMODEL%YRML_GCONF%YGFL%YCHEM(JLAYER-NSV_CHEMBEG_MSE+1)%RMOLMASS 
      ZSVM(KST:KEND,JLAYER) = PSVM(KST:KEND,JLAYER) * 1000.0 * XAVOGADRO * &
          PRHODREF(KST:KEND) / ZMOLMASS(JLAYER)
    END DO
  ELSE
    DO JLAYER=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
      !input unit assumed to be mol/mol (default case)
      ZSVM(KST:KEND,JLAYER) = PSVM(KST:KEND,JLAYER) * XAVOGADRO * PRHODREF(KST:KEND) / XMD
    END DO
  END IF

  DO JLAYER=NSV_AERBEG_MSE,NSV_AEREND_MSE ! MC : isn't aerosol unit kg/kg?
    ZSVM(KST:KEND,JLAYER) = PSVM(KST:KEND,JLAYER) * XAVOGADRO * PRHODREF(KST:KEND) / XMD
  END DO
  DO JLAYER=NSV_DSTBEG_MSE,NSV_DSTEND_MSE
    ZSVM(KST:KEND,JLAYER) = PSVM(KST:KEND,JLAYER) * XMOLARWEIGHT_DST*&
        PRHODREF(KST:KEND) / XMD
  END DO
  DO JLAYER=NSV_DSTDEPBEG_MSE,NSV_DSTDEPEND_MSE
    ZSVM(KST:KEND,JLAYER) = 0.
  END DO
ENDIF

!        1.14 Implicit coupling coeff.

IF(.NOT.LD_MPA) THEN
  ZCFBUV(KST:KEND) = SQRT( PCFBU(KST:KEND)**2 + PCFBV(KST:KEND)**2 )
  ZCFBQ (KST:KEND) = PCFBQ(KST:KEND) * ZRHOA(KST:KEND)
  ZCFAQ (KST:KEND) = PCFAQ(KST:KEND) * ZRHOA(KST:KEND)
ELSE
  ZCFBUV(KST:KEND) = 0._JPRB
  ZCFBQ (KST:KEND) = 0._JPRB
  ZCFAQ (KST:KEND) = PCFAQ(KST:KEND)
ENDIF


!*       2.     Call to surface monitor with 2D variables
!               -----------------------------------------

ILEN=KEND-KST+1

ZTSURF(KST:KEND) = 0._JPRB
ZZ0(KST:KEND)    = 0._JPRB
ZZ0H(KST:KEND)   = 0._JPRB
ZQSURF(KST:KEND) = 0._JPRB
ZIMPWET(KST:KEND,:) = 0._JPRB
ZIMPDRY(KST:KEND,:) = 0._JPRB

CALL COUPLING_SURF_ATM_N(YSURFEX_ARO_CUR, CPROGRAM_MSE, CCOUPLING, PTIMEC, &
   PTSTEP, IYEAR, IMONTH, IDAY, PRHGMT, ILEN, KSV, KSW,&
   ZTSUN(KST:KEND), ZZENITH(KST:KEND), ZZENITH2(KST:KEND), ZAZIM(KST:KEND),&
   ZZREF(KST:KEND), ZZREF(KST:KEND), PZS(KST:KEND),&
   ZU(KST:KEND), ZV(KST:KEND), ZQA(KST:KEND),&
   ZTA(KST:KEND), ZRHOA(KST:KEND),ZSVM(KST:KEND,:),&
   ZCO2(KST:KEND),ZIMPWET(KST:KEND,:),ZIMPDRY(KST:KEND,:), YSV,&
   ZRAIN(KST:KEND), ZSNOW(KST:KEND), PFRTHDS(KST:KEND),&
   PFRSODIDS(KST:KEND,:), PFRSOSCDS(KST:KEND,:), PSW_BANDS,&
   ZPS(KST:KEND), ZPA(KST:KEND),&
   ZSFTQ(KST:KEND) , ZSFTH(KST:KEND),&
   ZSFTS(KST:KEND,:),ZSFCO2(KST:KEND),&
   ZSFU(KST:KEND), ZSFV(KST:KEND),&
   PTSRAD(KST:KEND), PDIR_ALB(KST:KEND,:),&
   PSCA_ALB(KST:KEND,:), PEMIS(KST:KEND),&
   ZTSURF(KST:KEND), ZZ0(KST:KEND), ZZ0H(KST:KEND), ZQSURF(KST:KEND), &
   PCFAU(KST:KEND), ZCFBUV(KST:KEND),&
   PCFATH(KST:KEND), ZCFAQ(KST:KEND),&
   PCFBTH(KST:KEND), ZCFBQ(KST:KEND),&
   'OK')

!LWNETF(KST:KEND) = 0._JPRB

! Friction of components along slope axes (U: largest local slope axis, V: zero slope axis)

IF (LD_MPA) THEN

  WHERE (ZSFU(KST:KEND)/=XUNDEF .AND. ZWIND(KST:KEND)>0.)
    PSFU(KST:KEND) = - SQRT(ZSFU(KST:KEND)**2+ZSFV(KST:KEND)**2) *  &
       ZUA(KST:KEND) / ZWIND(KST:KEND) / PRHODREF(KST:KEND)
    PSFV(KST:KEND) = - SQRT(ZSFU(KST:KEND)**2+ZSFV(KST:KEND)**2) *  &
       ZVA(KST:KEND) / ZWIND(KST:KEND) / PRHODREF(KST:KEND)
  END WHERE

  !* conversion from H (W/m2) to w'Theta'

  WHERE (ZSFTH(KST:KEND)/=XUNDEF)
    PSFTH(KST:KEND) = ZSFTH(KST:KEND) / XCPD / PRHODREF(KST:KEND)
  END WHERE

  !* conversion from water flux (kg/m2/s) to w'rv'

  WHERE (ZSFTQ(KST:KEND)/=XUNDEF)
    PSFRV(KST:KEND) = ZSFTQ(KST:KEND) / PRHODREF(KST:KEND)
  ENDWHERE

ELSE

  IF (CCOUPLING == 'E') THEN
    WHERE (ZSFU(KST:KEND)/=XUNDEF .AND. ZWIND(KST:KEND)>0.)
      ZSFMOD(KST:KEND) = SQRT(ZSFU(KST:KEND)**2+ZSFV(KST:KEND)**2)
      PSFU(KST:KEND) = ZSFMOD(KST:KEND) * ZUA(KST:KEND) / ZWIND(KST:KEND)
      PSFV(KST:KEND) = ZSFMOD(KST:KEND) * ZVA(KST:KEND) / ZWIND(KST:KEND)
    END WHERE
  ELSEIF (CCOUPLING == 'I') THEN
    WHERE (ZSFU(KST:KEND)/=XUNDEF .AND. ZCFBUV(KST:KEND)>0.)
      ZSFMOD(KST:KEND) = SQRT(ZSFU(KST:KEND)**2+ZSFV(KST:KEND)**2)
      PSFU(KST:KEND) = ZSFMOD(KST:KEND) * PCFBU(KST:KEND) / ZCFBUV(KST:KEND)
      PSFV(KST:KEND) = ZSFMOD(KST:KEND) * PCFBV(KST:KEND) / ZCFBUV(KST:KEND)
    END WHERE
  ELSEIF (CCOUPLING == 'V') THEN
    WHERE (ZSFU(KST:KEND)/=XUNDEF .AND. ZWIND(KST:KEND)>0.)
      ZSFMOD(KST:KEND) = SQRT(ZSFU(KST:KEND)**2+ZSFV(KST:KEND)**2) / &
        ZWIND(KST:KEND)
      PSFU(KST:KEND) = ZSFMOD(KST:KEND) * PCFBU(KST:KEND) / &
        ( 1.0 - ZSFMOD(KST:KEND) * PCFAU(KST:KEND) )
      PSFV(KST:KEND) = ZSFMOD(KST:KEND) * PCFBV(KST:KEND) / &
        ( 1.0 - ZSFMOD(KST:KEND) * PCFAU(KST:KEND) )
    END WHERE
  ELSE
    CALL ABOR1('ARO_GROUND_PARAM: CCOUPLING WRONG VALUE!')
  ENDIF

  WHERE (ZSFTH(KST:KEND)/=XUNDEF)
    PSFTH(KST:KEND) = -ZSFTH(KST:KEND)
  ENDWHERE

  WHERE (ZSFTQ(KST:KEND)/=XUNDEF)
    PSFRV(KST:KEND) = -ZSFTQ(KST:KEND)
  ENDWHERE

ENDIF

IF (LELAM) THEN
  !  In E zone, pemis=0. so Tsrad -> +Inf. Initialisations for security for radiation
  ZDGUN = REAL(KDGUNG)
  ZDGUX = REAL(KDGUXG)
  ZDLUN = REAL(KDLUNG)
  ZDLUX = REAL(KDLUXG)

  DO JI=KST,KEND
    ZEMASK(JI)=MAX(0.0,SIGN(1.0,(PINDX(JI)-ZDLUN)*(ZDLUX-PINDX(JI))))&
     & *MAX(0.0,SIGN(1.0,(PINDY(JI)-ZDGUN)*(ZDGUX-PINDY(JI))))
  ENDDO

  DO  JI= KST, KEND
    IF (ZEMASK(JI)== 0.) THEN
      PTSRAD  (JI)=288.
      PEMIS   (JI)=0.99
      PDIR_ALB(JI,:)=0.1
      PSCA_ALB(JI,:)=0.1
    ENDIF
  ENDDO
ENDIF
!  conversion from scalar flux (kg/m2/s) to w'rsv'
IF(KSV .GT. 0.AND..NOT.LCHEM_ARPCLIM) THEN  !MM I don't understand in which general case these instructions are ok...
  DO JSV=1,KSV
    WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF)
      PSFSV(KST:KEND,JSV)=ZSFTS(KST:KEND,JSV)/PRHODREF(KST:KEND)
    ENDWHERE
  ENDDO
ENDIF

!  conversion unit of chemistry or aerosol flux 
IF ((LUSECHEM_MSE).AND.(LCH_SURFACE_FLUX_MSE)) THEN
!  (molec/m2/s) to (ppp.m.s-1)
  DO JSV=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
    WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF)
      PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * XMD / ( XAVOGADRO * PRHODREF(KST:KEND))
    END WHERE
    WRITE(20,*) 'FLUX ppp.m.s-1 : ',JSV-NSV_CHEMBEG_MSE+1,' =',&
        MINVAL(PSFSV(KST:KEND,JSV)),&
        MAXVAL(PSFSV(KST:KEND,JSV))
  END DO
ELSEIF (LLCHEMAERO_ARPC) THEN
  DO JSV=NSV_CHEMBEG_MSE,NSV_CHEMBEG_MSE+NAERO-1
    WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF) 
      ! (kg/m2/s) to (kg/m2/s) 
      PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV)
    END WHERE
  END DO
  IF (LCHEM_ARPCLIM.AND.TRIM(CHEM_SCHEME).EQ.'mocage') THEN
       DO JSV=NSV_CHEMBEG_MSE+NAERO,NSV_CHEMEND_MSE
          WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF) 
             ! (molec/m2/s) to (kg/m2/s) 
             PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * ZMOLMASS(JSV) / 1000.0 / XAVOGADRO
          END WHERE
       END DO
  ENDIF
  IF (LCHEM_ARPCLIM.AND.TRIM(CHEM_SCHEME).EQ.'arpclim_relacs') THEN
       ! Convert chemical surface fluxes from molecules_x m-2 s-1 to mol_x mol_air-1 kg_air m-2 s-1
       ! emissions fluxes converted in updclisfx into molecules_x m-2 s-1 to allow to add emission 
       ! and dry deposition fluxes 
       ! dry deposition velocity in m s-1 concentrations converted in molecules_x m-3 see above (ZSVM) 
       ! in aplpar chemical concentrations in mol mol-1 (=mol_x mol_air-1)
       DO JSV=NSV_CHEMBEG_MSE+NAERO,NSV_CHEMEND_MSE
          WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF)
             PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * & 
                  & (XMD/XAVOGADRO) 
          END WHERE
       END DO
  ENDIF
ELSEIF (LLCHEM_ARPC) THEN
  IF (LCHEM_ARPCLIM.AND.TRIM(CHEM_SCHEME).EQ.'mocage') THEN
       DO JSV=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
          WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF) 
             ! (molec/m2/s) to (kg/m2/s) 
             PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * ZMOLMASS(JSV) / 1000.0 / XAVOGADRO
          END WHERE
       END DO
  ENDIF
  IF (LCHEM_ARPCLIM.AND.TRIM(CHEM_SCHEME).EQ.'arpclim_relacs') THEN
       ! Convert chemical surface fluxes from molecules_x m-2 s-1 to mol_x mol_air-1 kg_air m-2 s-1
       ! emissions fluxes converted in updclisfx into molecules_x m-2 s-1 to allow to add emission 
       ! and dry deposition fluxes 
       ! dry deposition velocity in m s-1 concentrations converted in molecules_x m-3 see above (ZSVM) 
       ! in aplpar chemical concentrations in mol mol-1 (=mol_x mol_air-1)
       DO JSV=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
          WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF)
             PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * & 
                  & (XMD/XAVOGADRO) 
          END WHERE
       END DO
  ENDIF
ELSEIF (LLAERO_ARPC) THEN
  DO JSV=NSV_CHEMBEG_MSE,NSV_CHEMBEG_MSE+NAERO-1
    WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF) 
      ! (molec/m2/s) to (kg/m2/s) 
      PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV)
    END WHERE
  END DO
ELSE
  PSFSV(KST:KEND,NSV_CHEMBEG_MSE:NSV_CHEMEND_MSE) = 0.
END IF

!* conversion from aerosol flux (not tactic scheme case) (molec/m2/s) to (ppp.m.s-1)
IF (LORILAM_MSE) THEN
  DO JSV=NSV_AERBEG_MSE,NSV_AEREND_MSE
    WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF)
    PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * XMD / ( XAVOGADRO * PRHODREF(KST:KEND))
    END WHERE
  END DO
ELSE
    PSFSV(KST:KEND,NSV_AERBEG_MSE:NSV_AEREND_MSE) = 0.
END IF

!* conversion from dust flux (kg/m2/s) to (ppp.m.s-1)
IF (LDUST_MSE) THEN
  DO JSV=NSV_DSTBEG_MSE,NSV_DSTEND_MSE
    WHERE (ZSFTS(KST:KEND,JSV)/=XUNDEF)
      PSFSV(KST:KEND,JSV) = ZSFTS(KST:KEND,JSV) * XMD &
            & / (XMOLARWEIGHT_DST * PRHODREF(KST:KEND))
    END WHERE
  END DO
ELSE
    PSFSV(KST:KEND,NSV_DSTBEG_MSE:NSV_DSTEND_MSE) = 0.
END IF
IF (LDEPOS_DST_MSE(1)) PSFSV(KST:KEND,NSV_DSTDEPBEG_MSE:NSV_DSTDEPEND_MSE) = 0.


IF (LD_MPA) THEN
  IF (NSV_CO2_MSE /= 0) THEN
  !* conversion from CO2 flux (kg/m2/s) to ppm.m.s-1
    PSFCO2(KST:KEND) = ZSFCO2(KST:KEND) / &
                     (PRHODREF(KST:KEND) * 1E-6 * 44./28.8)
  ELSE
  !* conversion from CO2 flux (kg/m2/s) to  w'CO2'
    PSFCO2(KST:KEND) = ZSFCO2(KST:KEND) / PRHODREF(KST:KEND)
  END IF

ELSE
   PSFCO2(KST:KEND) = ZSFCO2(KST:KEND)
ENDIF
DEALLOCATE(YSV)

IF (LDXFUMSE) THEN
  IF (YSURFEX_ARO_CUR%U%CWATER == 'FLAKE' ) CALL CORRECT_TIME_FLAKE(YSURFEX_ARO_CUR, PTSTEP)
  IF (YSURFEX_ARO_CUR%U%CWATER == 'WATFLX') CALL CORRECT_TIME_WAT(YSURFEX_ARO_CUR, PTSTEP)
  IF (YSURFEX_ARO_CUR%U%CNATURE == 'ISBA') CALL CORRECT_TIME_ISBA(YSURFEX_ARO_CUR, PTSTEP)
  IF (YSURFEX_ARO_CUR%U%CSEA  == 'SEAFLX') CALL CORRECT_TIME_SEA(YSURFEX_ARO_CUR, PTSTEP)
  IF (YSURFEX_ARO_CUR%U%CTOWN == 'TEB') CALL CORRECT_TIME_TEB(YSURFEX_ARO_CUR, PTSTEP)
  CALL CORRECT_TIME_SURF(YSURFEX_ARO_CUR, PTSTEP)
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ARO_GROUND_PARAM',1,ZHOOK_HANDLE)

END SUBROUTINE ARO_GROUND_PARAM



