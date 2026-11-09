!OPTIONS XOPT(NOEVAL)
SUBROUTINE RECMWF(YDDIMV,YDMODEL,                      &
 & KIDIA   , KFDIA   , KLON     , KLEV     , KSW     , &
 & KOZOCL  , KMACCAERO, KAERO, &
 & PALBD   , PALBP   , PAPRS    , PAPRSF   ,           &
 & PNEB    , PQO3    , PQCO2    , PQCH4    , PQN2O   , &
 & PQNO2   , PQC11   , PQC12    , PQC22    , PQCL4   , &
 & PAER    , PAERO   , PAEROCPL , PDELP    , PEMIS    , &
 & PMU0M   , PQ      , PQSAT    , PQICE    , PQLIQ   , &
 & PS      , PRR     , PLSM     , PT       , PTS     , &
 & PEMTD   , PEMTU   , PTRSO    , PLWFC    , PLWFT   , &
 & PSWFC   , PSWFT   , PSFSWDIR , PSFSWDIF , PFSDNN  , &
 & PFSDNV  , PCTRSO  , PCEMTR   , PTRSOD   , PTRSODIR, &
 & PTRSODIF,PPIZA_DST,PCGA_DST  ,PTAUREL_DST,&
 & PAERINDS, PGELAM, PGEMU,      PSWDIR , PSWDIF  , &
 & PMU0LU , PALB, PVRMOON, &
 & PCLDROP, YSPP_RSWINHF, YSPP_RLWINHF)

!**** *RECMWF* - METEO-FRANCE RADIATION INTERFACE TO ECMWF RADIATION SCHEME

!**** PURPOSE :
!     ---------
!           SIMPLE INTERFACE TO RADLSW - RADLSWR - RADIATION_SCHEME

!**   INTERFACE.
!     ----------

!     EXPLICIT ARGUMENTS :
!        --------------------
! KIDIA    : START INDEX OF DATA IN KLON-LONG VECTOR
! KFDIA    : END   INDEX OF DATA IN KLON-LONG VECTOR
! KLON     : VECTOR LENGTH
! KLEV     : NUMBER OF LEVELS
! KOZOCL   : NOZOCL (Ozone clim choice)
! KMACCAERO: MACC AEROSOLS (clims or prognostics)
! KAERO    : NUMBER OF MACC AEROSOLS 
! PAER     : (KLON,KLEV ,6)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PAERO    : (KLON,KLEV ,KAERO)    ; AEROSOLS MACC MMR 
!                                     (from netcdf clim file or GFL)
! PAEROCPL : (KLON,KLEV ,16,6)  ; OPTICAL AEROSOL PROPERTIES FOR RADIATION 
! PAERINDS : (KLON,KLEV)        ; OPTICAL THICKNESS OF THE SULFATE AEROSOL (indirect effect only)
! PALBD    : (KLON,NSW)         ; DIFFUSE ALBEDO IN THE 2 SW INTERVALS
! PALBP    : (KLON,NSW)         ; PARALLEL ALBEDO IN THE 2 SW INTERVALS
! PAPRS    : (KLON,KLEV+1)      ; HALF LEVEL PRESSURE
! PAPRSF   : (KLON,KLEV )       ; FULL LEVEL PRESSURE
! PNEB     : (KLON,KLEV )       ; CLOUD FRACTIONAL COVER
! PQO3     : (KLON,0:KLEV )     ; OZONE MIXING RATIO (MASS)
! PQCO2    : (KLON,KLEV)        ; CO2 from MACC clims 
! PQCH4    : (KLON,KLEV)        ; from MACC clims 
! PQN2O    : (KLON,KLEV)        ; from MACC clims   
! PQNO2    : (KLON,KLEV)        ; from MACC clims
! PQC11    : (KLON,KLEV)        ; from MACC clims   
! PQC12    : (KLON,KLEV)        ; from MACC clims    
! PQC22    : (KLON,KLEV)        ; from MACC clims    
! PQCL4    : (KLON,KLEV)        ; from MACC clims
! PDELP    : (KLON,KLEV)        ; LAYER PRESSURE THICKNESS
! PEMIS    : (KLON)             ; SURFACE EMISSIVITY
! PMU0M    : (KLON)             ; MEAN SOLAR ANGLE
! PQ       : (KLON,KLEV )       ; SPECIFIC HUMIDITY PA/PA
! PQSAT    : (KLON,KLEV )       ; SATURATION SPECIFIC HUMIDITY PA/PA
! PQICE    : (KLON,KLEV )       ; ICE    WATER KG/KG
! PQLIQ    : (KLON,KLEV )       ; LIQUID WATER KG/KG
! PLSM     : (KLON)             ; LAND-SEA MASK
! PT       : (KLON,KLEV)        ; FULL LEVEL TEMPERATURE
! PTS      : (KLON)             ; SURFACE TEMPERATURE
! PPIZA_DST: (KLON,KLEV,NSW); Single scattering albedo of dust 
! PCGA_DST : (KLON,KLEV,NSW); Assymetry factor for dust 
! PTAUREL_DST: (KLON,KLEV,NSW); Optical depth of dust relative to at 550nm

!     ==== OUTPUTS ===
! PEMTD (KLON,KLEV+1)         ; TOTAL DOWNWARD LONGWAVE EMISSIVITY
! PEMTU (KLON,KLEV+1)         ; TOTAL UPWARD   LONGWAVE EMISSIVITY
! PTRSO (KLON,KLEV+1)         ; TOTAL SHORTWAVE TRANSMISSIVITY
! PCTRSO(KLON,2)              ; CLEAR-SKY SHORTWAVE TRANSMISSIVITY
! PCEMTR(KLON,2)              ; CLEAR-SKY NET LONGWAVE EMISSIVITY
! PTRSOD(KLON)                ; TOTAL-SKY SURFACE SW TRANSMISSITY
! PLWFC (KLON,2)              ; CLEAR-SKY LONGWAVE FLUXES
! PLWFT (KLON,KLEV+1)         ; TOTAL-SKY LONGWAVE FLUXES
! PSWFC (KLON,2)              ; CLEAR-SKY SHORTWAVE FLUXES
! PSWFT (KLON,KLEV+1)         ; TOTAL-SKY SHORTWAVE FLUXES
! PTRSODIR(KLON,NSW)          ; DIRECT SURFACE SW TRANSMISSITY
! PTRSODIF(KLON,NSW)          : DIFFUSE SURFACE SW TRANSMISSITY

!   === OPTIONAL ===
! PCLDROP(KLON,KLEV)          ; Cloud Droplet number conc. from n.r.t. aerosols

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!     SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHORS.
!     --------
!     ORIGINAL BY  B. RITTER   *ECMWF*        83-10-13
!     REWRITING FOR IFS BY J.-J. MORCRETTE    94-11-15

!     MODIFICATIONS.
!     --------------
!     2003-09: Y.   Bouteloup              : Duplication of RFMR to use present (cy25) ECMWF radiation scheme
!     2004-09: F.   Bouyssel               : Use of 6 aerosols & introduce NSW
!     2018-04: Y.   Seity                  : 4 New arguments for AROME 
!     2005-10: Y.   Seity                  : 3 optional arguments for dust optical properties 
!     2006-07: J.J. Morcrette              : PP of clear-sky PAR and TOA incident solar radiation (ECMWF)
!     2011-02: A.   Voldoire               : add PAERINDS for sulfate indirect effect computation
!     2011-06: M.   Jerczynski             : some cleaning to meet norms
!     2014-05: J.M. Piriou                 : force qv used by radiation to be larger to some factor of qsat.
!     2014-11: M.Benamara  + Y. Bouteloup  : Call to RADLSWR if LSRTM
!     2017-09: M.Bouzghaia + Y.  Bouteloup : Call to ECRAD
!     2018-09: M.Bouzghaia + Y.  Bouteloup : ECRAD adaptation to MF + Cleanings
!     2019-10: Y. Bouteloup : Phasing to cy47t0
!     2020-12: U. Andrae : Introduce SPP for HARMONIE-AROME
!     2021-01: D. Martin Add cloud droplet number conc. from nrt aerosols
!     2022-01: Y. Bouteloup and Yann Seity : Fixes to run Arome with EcRad and cleaning of commented lines
!     21-Jun-2022 R. El Khatib portability fix
!     2023-01, M. Michou : modifications so that chemistry interacts with ecRad
!     2023-09, P.Nabat : add PAEROCPL for TACTIC aerosol-radiation coupling
!-----------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE YOMDIMV    , ONLY : TDIMV
USE PARKIND1   , ONLY : JPIM     ,JPRB
USE YOMCST     , ONLY : RMCH4  ! RMO3 , RMCO2 , RMCH4 , RMCO, RMHCHO , RMSO2, RMSO4
USE YOMHOOK    , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE SPP_MOD_TYPE, ONLY : TSPP_CONFIG_TYPE,APPLY_SPP 

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KOZOCL
INTEGER(KIND=JPIM),INTENT(IN)    :: KMACCAERO
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEB(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQO3(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCO2(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCH4(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQN2O(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQNO2(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQC11(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQC12(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQC22(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCL4(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,KLEV,6) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV,KAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAEROCPL(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NTLW,6)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERINDS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0M(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQICE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLIQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMTD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMTU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRSO(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTRSO(KLON,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCEMTR(KLON,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRSOD(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWFC(KLON,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWFT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWFC(KLON,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWFT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSFSWDIR(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSFSWDIF(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDNN(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSDNV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRSODIR(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRSODIF(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PS(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSWDIR(KLON,KSW)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSWDIF(KLON,KSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0LU(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PALB(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVRMOON(KLON)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PCLDROP(KLON,KLEV)
TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT), OPTIONAL :: YSPP_RSWINHF, YSPP_RLWINHF 

!*       0.2   LOCAL ARRAYS.
!              -------------
REAL(KIND=JPRB) :: ZTH(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZRAER  (KLON,6,KLEV)
REAL(KIND=JPRB) :: ZRCLC  (KLON,KLEV)
REAL(KIND=JPRB) :: ZRMU0  (KLON)
REAL(KIND=JPRB) :: ZQLWP  (KLON,KLEV ) , ZQIWP (KLON,KLEV )
REAL(KIND=JPRB) :: ZPQO3  (KLON,KLEV)
REAL(KIND=JPRB) :: ZQS    (KLON,KLEV)
REAL(KIND=JPRB) :: ZQ     (KLON,KLEV)
REAL(KIND=JPRB) :: ZCLDROP(KLON,KLEV)
REAL(KIND=JPRB) :: ZNBAS  (KLON) , ZNTOP  (KLON)
REAL(KIND=JPRB) :: ZCCNL  (KLON) , ZCCNO  (KLON)
REAL(KIND=JPRB) :: ZPFRSODC(KLON)
REAL(KIND=JPRB) :: ZPFRTEDC(KLON)
INTEGER(KIND=JPIM) :: IAERO    
REAL(KIND=JPRB) :: ZRCH4(KLON,KLEV), ZRCO2(KLON,KLEV), ZRN2O(KLON,KLEV), ZRNO2(KLON,KLEV)
REAL(KIND=JPRB) :: ZRC11(KLON,KLEV), ZRC12(KLON,KLEV), ZRC22(KLON,KLEV), ZRCL4(KLON,KLEV)
REAL(KIND=JPRB) ::  ZPAERO  (KLON,KLEV,KAERO)
REAL(KIND=JPRB) ::  ZPFDIR  (KLON) , ZPFDIF   (KLON), ZPCDIR         (KLON)
REAL(KIND=JPRB) ::  ZFLUX_UV(KLON) , ZFLUX_PAR(KLON), ZFLUX_PAR_CLEAR(KLON) 
!  output of radlswr

REAL(KIND=JPRB) :: ZFRTED(KLON)
REAL(KIND=JPRB) :: ZLWDERIVATIVE(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZSWDIFFUSEBAND(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZSWDIRECTBAND(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)

!  output of radlsw

REAL(KIND=JPRB) :: ZEMIT  (KLON)
REAL(KIND=JPRB) :: ZFCT   (KLON,KLEV+1)
REAL(KIND=JPRB) :: ZFCS   (KLON,KLEV+1)
REAL(KIND=JPRB) :: ZFRSOD (KLON),ZSUDU(KLON)
REAL(KIND=JPRB) :: ZPARF  (KLON),ZUVDF(KLON),ZPARCF(KLON),ZTINCF(KLON)
REAL(KIND=JPRB) :: ZRSWINHF(KLON),ZRLWINHF(KLON)

INTEGER(KIND=JPIM) :: JK, JL, JSW
REAL(KIND=JPRB)    :: ZCRAE, ZEMIW(KLON),ZFRAC, ZSW
REAL(KIND=JPRB)    :: ZSPECTRALEMISS(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NLWEMISS)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

#include "radlsw.intfb.h"
#include "radlswr.intfb.h"
#include "radiation_scheme.intfb.h"
#include "acralu.intfb.h"
! ================================================================================


IF (LHOOK) CALL DR_HOOK('RECMWF',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & YDARPHY => YDMODEL%YRML_PHY_MF%YRARPHY  , YDERAD  => YDMODEL%YRML_PHY_RAD%YRERAD  , &
 & YDPHY3  => YDMODEL%YRML_PHY_MF%YRPHY3   , YDERDI  => YDMODEL%YRML_PHY_RAD%YRERDI  , &
 & YDPHY   => YDMODEL%YRML_PHY_MF%YRPHY    , YDRIP   => YDMODEL%YRML_GCONF%YRRIP     , & 
 & YDTOPH  => YDMODEL%YRML_PHY_MF%YRTOPH   , YDEAERD => YDMODEL%YRML_PHY_RAD%YREAERD , &
 & YDPARAR => YDMODEL%YRML_PHY_MF%YRPARAR  , YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM)
 
ASSOCIATE( LUSEPRE2017RAD=>YDERAD%LUSEPRE2017RAD      , &
 & LSRTM   =>YDERAD%LSRTM   , LRRTM  =>YDERAD%LRRTM   , &
 & NAER    =>YDERAD%NAER    , NSW    =>YDERAD%NSW     , &
 & REPH2O  =>YDERDI%REPH2O  , RRAE   =>YDERDI%RRAE    , &
 & RSWINHF=>YDERAD%RSWINHF, RLWINHF=>YDERAD%RLWINHF   , &
 & LRDUST  =>YDARPHY%LRDUST , GQVTOP =>YDPARAR%GQVTOP , &
 & XCQVR   =>YDPARAR%XCQVR  , RII0   =>YDPHY3%RII0    , &
 & LAERONRTR=>YDEAERATM%LAERONRTR)

!*       1.    PREPARATORY WORK
!              ----------------

!*       1.1    LOCAL CONSTANTS
!                ---------------
IF (KOZOCL==4.OR.KOZOCL==5.OR.KOZOCL==6) THEN ! note that H2O or PQ is not in the list below
! coupling of radiative and chemical water vapor is under CHEM_LCHQ key  
  ZRCH4(:,:)   = PQCH4(:,:)
  ZRN2O(:,:)   = PQN2O(:,:)
  ZRNO2(:,:)   = PQNO2(:,:)
  ZRC11(:,:)   = PQC11(:,:)
  ZRC12(:,:)   = PQC12(:,:)
  ZRC22(:,:)   = PQC22(:,:)
  ZRCL4(:,:)   = PQCL4(:,:)
  ZRCO2(:,:)   = PQCO2(:,:)
ELSE
  ZRCH4(:,:)   = YDERDI%RCH4
  ZRN2O(:,:)   = YDERDI%RN2O
  ZRNO2(:,:)   = YDERDI%RNO2
  ZRC11(:,:)   = YDERDI%RCFC11
  ZRC12(:,:)   = YDERDI%RCFC12
  ZRC22(:,:)   = YDERDI%RCFC22
  ZRCL4(:,:)   = YDERDI%RCCL4
  ZRCO2(:,:)   = YDERDI%RCARDI
ENDIF
! IF (KAERO>0.OR.(TRIM(YDMODEL%YRML_CHEM%YRCOMPO%AERO_SCHEME)=="tactic".AND.YDMODEL%YRML_PHY_MF%YRPHY%LACLIB).OR.YDMODEL%YRML_PHY_RAD%YREAERATM%LAERONRTR) THEN
IF (KAERO>0 .OR. LAERONRTR) THEN
  ZPAERO(:,:,:)= PAERO(:,:,:)
ELSE
  ZPAERO(:,:,:)= 0.0_JPRB
ENDIF
ZRAER(:,:,:) = 0.0_JPRB
IAERO        = KAERO
ZNBAS(:)     = 1.0_JPRB
ZNTOP(:)     = 1.0_JPRB
ZCCNL(:)     = YDERAD%RCCNLND
ZCCNO(:)     = YDERAD%RCCNSEA
ZCRAE        = RRAE*(RRAE+2.0_JPRB)
IF (PRESENT(PCLDROP)) THEN
  ZCLDROP(:,:) = PCLDROP(:,:)
ENDIF

  !*       2.1    FULL-LEVEL QUANTITIES
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZPQO3(JL,JK)=PQO3(JL,JK)*PDELP(JL,JK)
    ZRCLC(JL,JK)=MAX( 0.0_JPRB ,MIN( 1.0_JPRB ,PNEB(JL,JK)))
    IF(YDPARAR%LICERAD)THEN ! Use fraction of ice/snow/graupel near one instead of
!      fraction as cloud cover
       ZFRAC=PQICE(JL,JK)/(PQLIQ(JL,JK) + PQICE(JL,JK) + 1.0E-15_JPRB)
       ZRCLC(JL,JK) = ZFRAC + (1.0_JPRB-ZFRAC)*ZRCLC(JL,JK) ! Higher cloud condensate
!      fraction when more cloudice/snow
    ENDIF
    IF (ZRCLC(JL,JK) > YDERDI%REPCLC) THEN
      ZQLWP(JL,JK)=PQLIQ(JL,JK)
      ZQIWP(JL,JK)=PQICE(JL,JK)
    ELSE
      ZQLWP(JL,JK)=REPH2O*ZRCLC(JL,JK)
      ZQIWP(JL,JK)=REPH2O*ZRCLC(JL,JK)
    ENDIF
    ZQS (JL,JK)=MAX(2.0_JPRB*REPH2O,PQSAT(JL,JK))
    ZQ  (JL,JK)=MAX(REPH2O,MIN(PQ(JL,JK),ZQS(JL,JK)*(1.0_JPRB-REPH2O)))
    IF(YDPARAR%LQVTOP) THEN
      ! qv in input to radiation is set at least to XCQVR*qsat, in order to have
      ! greater qv between 100 and 200 hPa, and thus larger radiative heating at
      ! those levels. At the top of atmosphere, the value GQVTOP is used.
      ! The transition between GQVTOP and  XCQVR*qsat is done by ZFRAC,
      ! depending on pressure.
      ZFRAC=MIN(1._JPRB,PAPRSF(JL,JK)/YDPARAR%GQVPLIM)
      ZQ(JL,JK)=ZFRAC*MAX(XCQVR*PQSAT(JL,JK),ZQ(JL,JK))+(1._JPRB-ZFRAC)*GQVTOP
    ENDIF
  ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
        ZSPECTRALEMISS(JL,:)=PEMIS(JL)
        ZEMIW(JL)=PEMIS(JL)
  ENDDO

  IF (NAER == 0) THEN
    ZRAER=YDEAERD%RCAEROS
  ELSE
    DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      if (YDPHY%LAEROLAN) ZRAER(JL,1,JK)=PAER(JL,JK,1)
      if (YDPHY%LAEROSEA) ZRAER(JL,2,JK)=PAER(JL,JK,2)
      if (YDPHY%LAERODES) ZRAER(JL,3,JK)=PAER(JL,JK,3)
      if (YDPHY%LAEROSOO) ZRAER(JL,4,JK)=PAER(JL,JK,4)
      if (YDPHY%LAEROVOL) ZRAER(JL,5,JK)=YDEAERD%RCAEROS
      if (YDPHY%LAEROSUL) ZRAER(JL,6,JK)=PAER(JL,JK,6)
    ENDDO
    ENDDO
  ENDIF

  !*       2.2    HALF-LEVEL TEMPERATURE

  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      ZTH(JL,JK)=&
       & (PT(JL,JK-1)*PAPRSF(JL,JK-1)*(PAPRSF(JL,JK)-PAPRS(JL,JK))&
       & +PT(JL,JK)*PAPRSF(JL,JK)*(PAPRS(JL,JK)-PAPRSF(JL,JK-1)))&
       & *(1.0_JPRB/(PAPRS(JL,JK)*(PAPRSF(JL,JK)-PAPRSF(JL,JK-1))))  
    ENDDO
  ENDDO  
  DO JL=KIDIA,KFDIA
    ZTH(JL,KLEV+1)=PTS(JL)
    ZTH(JL,1)=PT(JL,1)-PAPRSF(JL,1)*(PT(JL,1)-ZTH(JL,2))&
    & /(PAPRSF(JL,1)-PAPRS(JL,2))  
  ENDDO

  !*       3.1     SOLAR ZENITH ANGLE IS EARTH'S CURVATURE CORRECTED

  DO JL=KIDIA,KFDIA
    IF (PMU0M(JL) > 1.E-10_JPRB) THEN
      ZRMU0(JL)=RRAE/(SQRT(PMU0M(JL)**2+ZCRAE)-PMU0M(JL))
    ELSE
      ZRMU0(JL)= RRAE/SQRT(ZCRAE)
    ENDIF
  ENDDO

! Setup SPP patterns
  IF (PRESENT(YSPP_RSWINHF)) THEN
    IF (YSPP_RSWINHF%LPERT) THEN
      CALL APPLY_SPP(YSPP_RSWINHF,KLON,KIDIA,KFDIA,RSWINHF,ZRSWINHF)
    ELSE
      ZRSWINHF(KIDIA:KFDIA) = RSWINHF
    ENDIF
  ELSE
    ZRSWINHF(KIDIA:KFDIA) = RSWINHF
  ENDIF
  IF (PRESENT(YSPP_RLWINHF)) THEN
    IF (YSPP_RLWINHF%LPERT) THEN
      CALL APPLY_SPP(YSPP_RLWINHF,KLON,KIDIA,KFDIA,RLWINHF,ZRLWINHF)
    ELSE
      ZRLWINHF(KIDIA:KFDIA) = RLWINHF
    ENDIF
  ELSE
    ZRLWINHF(KIDIA:KFDIA) = RLWINHF
  ENDIF

  !*       2.    CALL TO THE RADIATION SCHEME
  !              ----------------------------

  IF (LSRTM .and. LUSEPRE2017RAD) THEN
  !**** 2.a/ CALL TO RADIATION SCHEME WITH RRTM LW+SW  (SRTM)

    CALL RADLSWR(                                                              &
  ! ---- INPUTS ----------------------------------------------------------------
     & YDDIMV  , YDMODEL , KIDIA    , KFDIA   , KLON , KLEV   , IAERO  , RII0  ,        &
     & ZRAER    , ZPAERO  , PALBD   , PALBP  , PAPRS  , PAPRSF ,                &
     & ZCCNL   , ZCCNO   , PGELAM  , PGEMU  ,                                  &
     & ZRCO2   , ZRCH4   , ZRN2O   , ZRNO2  , ZRC11  , ZRC12  , ZRC22 , ZRCL4 ,&
     & PNEB   , PDELP     , PEMIS   , ZEMIW  , PLSM   , ZRMU0  , ZPQO3 ,        &
     & ZQ      , ZQIWP   , PS   , ZQLWP  , PRR  ,                         &
     & ZTH     , PT     , PTS     ,                                           &
  ! ---- OUTPUTS ---------------------------------------------------------------
     &  ZEMIT, ZFCT  , PLWFT , ZFCS , PSWFT  , ZFRSOD, &
     &  ZFRTED, ZPFRSODC, ZPFRTEDC ,&
     &  ZSUDU, ZUVDF , ZPARF, ZPARCF, ZTINCF, &
     &  ZPFDIR, ZPFDIF, ZPCDIR, &
     &  ZLWDERIVATIVE, ZSWDIFFUSEBAND, ZSWDIRECTBAND )

  ELSEIF (LRRTM .and. .not. LSRTM .and. LUSEPRE2017RAD) THEN
  !**** 2.b/ CALL TO RADIATION SCHEME WITH RRTM LW ONLY

    CALL RADLSW(                                                                 &
    ! ---- INPUTS ------------------------------------------------------------------
     & YDDIMV,YDMODEL%YRML_PHY_RAD,YDMODEL%YRML_PHY_SLIN%YREPHLI,YDMODEL%YRML_PHY_SLIN%YRPHNC,YDMODEL%YRML_PHY_EC%YRECLD, &
     & YDMODEL%YRML_PHY_MF%YRPARAR, &
     & KIDIA    , KFDIA     , KLON  , KLEV    , YDERAD%NMODE  , NAER ,                  &
     & RII0    ,                                                                 &
     & ZRAER    , PALBD    , PALBP   , PAPRS   , PAPRSF ,                         &
     & ZCCNL   , ZCCNO    ,                                                      &
     & YDERDI%RCARDI , ZRCLC    , PDELP     , PEMIS   , ZEMIW  , PLSM  , ZRMU0 , ZPQO3,  &
     & ZQ      , ZQIWP    , ZQLWP   , ZQS     ,                                  &
     & ZTH     , PT      , PTS     , ZNBAS   , ZNTOP  ,                         &
    ! ---- OUTPUTS -----------------------------------------------------------------
     & ZEMIT , ZFCT   , PLWFT    , ZFCS    , PSWFT,&
     & ZFRSOD, ZSUDU  , ZPARF   , ZUVDF   , ZPARCF, ZTINCF, PSFSWDIR,&
     & PSFSWDIF,PFSDNN,PFSDNV   ,&
     & LRDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST,PAERINDS, &
     & ZRSWINHF,ZRLWINHF,ZCLDROP)

  ELSEIF ( .not. LUSEPRE2017RAD ) THEN
  !**** 2.c/ CALL TO MODULAR RADIATION SCHEME ECRAD

    CALL RADIATION_SCHEME(YDMODEL                                 , &
    ! ---- INPUTS -------------------------------------------------------
     & KIDIA     , KFDIA     , KLON     , KLEV       , IAERO      , &
     & RII0     , ZRMU0    , PTS        , PALBD      , PALBP      , &
     & ZSPECTRALEMISS      ,                                        &
     & ZCCNL    , ZCCNO    ,                                        &
     & PGELAM   , PGEMU    , PLSM       ,                           &
     & PAPRSF   , PT       ,                                        &
     & PAPRS    , ZTH      ,                                        &
     & ZQ       , ZRCO2    , ZRCH4      , ZRN2O      , ZRNO2      , &
     & ZRC11    , ZRC12    , ZRC22      , ZRCL4      , ZPQO3      , &
     & ZRCLC    , ZQLWP    , ZQIWP      , PRR        , PS         , &
     & ZRAER    , ZPAERO   , PAEROCPL   ,                           &
   ! ---- OUTPUTS ------------------------------------------------------
     & PSWFT    , PLWFT    , ZFCS       , ZFCT       ,              &
     & ZFRSOD   , ZFRTED   , ZPFRSODC   , ZPFRTEDC   ,              &
     & ZPFDIR   , ZPCDIR   , ZSUDU      ,                           &
     & ZFLUX_UV , ZFLUX_PAR, ZFLUX_PAR_CLEAR         ,              &
     & ZTINCF   , ZEMIT    ,                                        &
     & ZLWDERIVATIVE       ,                                        &
     & ZSWDIFFUSEBAND      ,  ZSWDIRECTBAND )

    ZPFDIF(KIDIA:KFDIA) = ZFRSOD(KIDIA:KFDIA) - ZPFDIR(KIDIA:KFDIA)
  ENDIF

  DO JL=KIDIA,KFDIA
    PSWFC (JL, 1)=ZFCS(JL,     1)
    PSWFC (JL, 2)=ZFCS(JL,KLEV+1)
    PLWFC (JL, 1)=ZFCT(JL,     1)
    PLWFC (JL, 2)=ZFCT(JL,KLEV+1)
  ENDDO

  PEMTD=PLWFT
  PEMTU=PLWFT
  PCEMTR=PLWFC


  IF (LSRTM .AND. .NOT. LUSEPRE2017RAD) THEN   ! EcRad case
    DO JSW=1,NSW
      DO JL=KIDIA,KFDIA
        PTRSODIR(JL,JSW)=ZSWDIRECTBAND(JL,JSW)*(1.-PALBP(JL,JSW))/(RII0*ZRMU0(JL))
        PTRSODIF(JL,JSW)=ZSWDIFFUSEBAND(JL,JSW)*(1.-PALBD(JL,JSW))/(RII0*ZRMU0(JL))
      ENDDO
    ENDDO      
  ELSEIF (LSRTM .AND. LUSEPRE2017RAD) THEN     ! LSRTM / RRTM case
    ZSW = REAL(NSW,JPRB)
    DO JSW=1,NSW
       DO JL=KIDIA,KFDIA
         PTRSODIR(JL,JSW)=ZPFDIR(JL)*(1.-PALBP(JL,JSW))/ZSW/(RII0*ZRMU0(JL))
         PTRSODIF(JL,JSW)=ZPFDIF(JL)*(1.-PALBD(JL,JSW))/ZSW/(RII0*ZRMU0(JL))
      ENDDO
    ENDDO      
  ELSE                                         ! Other cases (SW and RRTM)
    DO JSW=1,NSW
       DO JL=KIDIA,KFDIA
         PTRSODIR(JL,JSW)=PSFSWDIR(JL,JSW)*(1.-PALBP(JL,JSW))/(RII0*ZRMU0(JL))
         PTRSODIF(JL,JSW)=PSFSWDIF(JL,JSW)*(1.-PALBD(JL,JSW))/(RII0*ZRMU0(JL))
       ENDDO
    ENDDO
  ENDIF

!*         4.2     TRANSFORM FLUXES TO MODEL HISTORICAL VARIABLES

  DO JK=1,KLEV+1
    DO JL=KIDIA,KFDIA
      PTRSO(JL,JK)=PSWFT(JL,JK)/(RII0*ZRMU0(JL))
    ENDDO
  ENDDO

  DO JK=1,2
    DO JL=KIDIA,KFDIA
      PCTRSO(JL,JK)=PSWFC(JL,JK)/(RII0*ZRMU0(JL))
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    PTRSOD(JL)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZFRSOD(JL)/(RII0*ZRMU0(JL))))
  ENDDO

  !*         7.3   RECONSTRUCT FLUXES FOR DIAGNOSTICS

  DO JL=KIDIA,KFDIA
    IF (PMU0M(JL) < 1.E-10_JPRB) THEN
      ZRMU0(JL)=0.0_JPRB
      PSFSWDIR(JL,:)=0.0_JPRB
      PSFSWDIF(JL,:)=0.0_JPRB
    ENDIF
  ENDDO

  IF (YDARPHY%LMSE) THEN
    IF (NSW==6) THEN
      DO JSW=1,NSW
        PSWDIR(KIDIA:KFDIA,JSW)= PTRSODIR(KIDIA:KFDIA,JSW)
        PSWDIF(KIDIA:KFDIA,JSW)= PTRSODIF(KIDIA:KFDIA,JSW)
      ENDDO
    ELSEIF (NSW==1) THEN
      PSWDIR(KIDIA:KFDIA,1)=PFSDNN(KIDIA:KFDIA)+PFSDNV(KIDIA:KFDIA)
      PSWDIF(KIDIA:KFDIA,1)=0.0_JPRB
    ENDIF
  ENDIF


  IF ( YDPHY%LRAYLU )  CALL ACRALU ( YDRIP,YDPHY3,KIDIA,KFDIA,KLON,YDTOPH%NTRADI,KLEV,&
                       & PAPRS, PAPRSF, PDELP, PNEB, PQ, ZRCO2, PQICE, PQLIQ, PQO3, PT,&
                       & PALB,PMU0LU, PVRMOON, PAER)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RECMWF',1,ZHOOK_HANDLE)
END SUBROUTINE RECMWF
