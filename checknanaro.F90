!OPTIONS XOPT(NOEVAL)
SUBROUTINE CHECKNANARO(YDRIP,YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,KLEV,KRR,KSTEP,&
  &PAPHI,PAPHIF,PAPRS,PAPRSF,PGELAM,PGEMU,PMU0,PLSM,PT,PQ,PQSAT,PTS,&
  &PTENDT, PTENDR, PTENDW, PTENDTKE,&
  &PFRTH, PFRSO, PFRTHDS, PFRSODS,&
  &PCLCH,PCLCL,PCLCM,PCLCT,PFPLSL,PFPLSN,PFPLSG,PSTRTU,PSTRTV,&
  &PFCS,PFCLL,PFCLN,PFEVL,PFEVN,ERRLOC)


!**** *CHECKNAN * - Check if fluxes are different from NaN.

!**   Interface.
!     ----------
!        *CALL* *CHECKNAN*

!-----------------------------------------------------------------------
! -   ARGUMENTS D'ENTREE.
! -   INPUT ARGUMENTS.
!     -------------------

! - DIMENSIONS.

! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
! KLON : HORIZONTAL DIMENSION                  (NPROMA IN *CPG*).
! KLEV : END OF VERTICAL LOOP AND VERTICAL DIMENSION(NFLEVG IN *CPG*).
! KRR  : NUMBER OF MICROPHYS. VARIABLES

! - 2D (0:KLEV) .

! PAPHI      : GEOPOTENTIAL ON HALF-LEVELS.
! PAPRS      : PRESSURE ON HALF-LEVELS.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIAL ON FULL LEVELS.
! PAPRSF     : PRESSURE ON FULL LEVELS.
! PT         : TEMPERATURE.
! PQ         : SPECIFIC HUMIDITY OF WATER VAPOUR.
! PQSAT      : SATURATION SPECIFIC HUMIDITY OF WATER VAPOUR.

! - 1D .
! PLSM       : LAND/SEA MASK.
! PTS        : SURFACE LAYER TEMPERATURE.
!-----------------------------------------------------------------------

! -   IMPLICIT ARGUMENTS.
!     ---------------------

!-----------------------------------------------------------------------

!     Externals.
!     ----------

!     Method.
!     -------

!     Author.
!     -------
!     2022-10-10, Y. Seity (From Checknan developped by J.-M. Piriou)

!     Modifications.
!     --------------
!     2023-11-08, R. H. Myhre Remove tendu and tendv that are longer in use
!                             (gone with LINTFLEX)

!     ------------------------------------------------------------------

USE IEEE_ARITHMETIC, ONLY : ISNAN => IEEE_IS_NAN
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMPHY0  , ONLY : TPHY0
USE YOMPHY2  , ONLY : TPHY2
USE YOMLUN_IFSAUX, ONLY : NULERR
USE YOMRIP0  , ONLY : NINDAT
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY0)       ,INTENT(IN) :: YDPHY0
TYPE(TPHY2)       ,INTENT(IN) :: YDPHY2
TYPE(TRIP)        ,INTENT(IN) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KRR
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQSAT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTS(KLON)

REAL(KIND=JPRB)   ,INTENT(IN), DIMENSION (KLON,KLEV,KRR) :: PTENDR
REAL(KIND=JPRB)   ,INTENT(IN), DIMENSION (KLON,KLEV) :: PTENDT, PTENDW, PTENDTKE
REAL(KIND=JPRB)   ,INTENT(IN), DIMENSION (KLON,0:KLEV) ::PFRSO,PFRTH,PFPLSL,&
                                         &PFPLSN,PFPLSG,PSTRTU,PSTRTV
REAL(KIND=JPRB)   ,INTENT(IN), DIMENSION (KLON) ::PFRTHDS,PFRSODS,PCLCH,PCLCL,&
                                 &PCLCM,PCLCT,PFCS,PFCLL,PFCLN,PFEVL,PFEVN

CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ERRLOC

REAL(KIND=JPRB), DIMENSION (KLON,0:KLEV) :: ZSUM

LOGICAL ::LLPTENDT, LLPTENDTKE, LLPTENDR,&
 &LLPFRSO,LLPFRTH,LLPFPLSL,&
 &LLPFPLSN,LLPFPLSG,LLPSTRTU,LLPSTRTV,LLPFRTHDS,LLPFRSODS,LLPCLCH,LLPCLCL,&
 &LLPCLCM,LLPCLCT,LLPFCS,LLPFCLL,LLPFCLN,LLPFEVL,LLPFEVN


REAL(KIND=JPRB) :: ZLOCST,ZCONRD,ZLON,ZLAT,ZHR
INTEGER(KIND=JPIM), PARAMETER :: JPERRN=100
INTEGER(KIND=JPIM) :: ILEV(JPERRN),ILON(JPERRN),JLEV,JLON,IERRN,JERRN,JRR
CHARACTER(LEN=200) :: CLVAR(JPERRN)
LOGICAL :: LLOK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHECKNANARO',0,ZHOOK_HANDLE)
ASSOCIATE(RSTATI=>YDRIP%RSTATI,RHGMT=>YDRIP%RHGMT)

!-------------------------------------------------
! Sum all fluxes: the sum vectorises (fast computation).
! The ISNAN test will be performed on the sum, and will performed on the
! individuals only if the sum contains NaN.
! This approach optimizes the computation time.
!-------------------------------------------------
DO JLEV=0,KLEV
  DO JLON=KIDIA,KFDIA
    ZSUM(JLON,JLEV)=0._JPRB
  ENDDO
ENDDO
DO JLEV=0,KLEV
  DO JLON=KIDIA,KFDIA
    ZSUM(JLON,JLEV)=ZSUM(JLON,JLEV)+ &
   & PFRSO(JLON,JLEV)+PFRTH(JLON,JLEV)+PFPLSL(JLON,JLEV)+&
   & PFPLSN(JLON,JLEV)+PFPLSG(JLON,JLEV)+&
   & PSTRTU(JLON,JLEV)+PSTRTV(JLON,JLEV)
  ENDDO
ENDDO
DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    ZSUM(JLON,JLEV)=ZSUM(JLON,JLEV)+ &
   &  PTENDT(JLON,JLEV) + PTENDW(JLON,JLEV) + PTENDTKE(JLON,JLEV)
  ENDDO
ENDDO
DO JRR=1,KRR
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZSUM(JLON,JLEV)=ZSUM(JLON,JLEV)+PTENDR(JLON,JLEV,JRR) 
    ENDDO
  ENDDO
ENDDO
DO JLON=KIDIA,KFDIA
  ZSUM(JLON,KLEV)=ZSUM(JLON,KLEV)+PFRTHDS(JLON)+PFRSODS(JLON)+&
  &PCLCH(JLON)+PCLCL(JLON)+PCLCM(JLON)+PCLCT(JLON)+&
  &PFCS(JLON)+PFCLL(JLON)+PFCLN(JLON)+PFEVL(JLON)+&
  &PFEVN(JLON)
ENDDO

LLOK=.TRUE.
DO JLEV=0,KLEV
  DO JLON=KIDIA,KFDIA
    IF(ISNAN(ZSUM(JLON,JLEV))) LLOK=.FALSE.
  ENDDO
ENDDO
IF(.NOT.LLOK) THEN
  !-------------------------------------------------
  ! Test fluxes individually.
  !-------------------------------------------------
  LLPTENDT=.TRUE.
  LLPTENDTKE=.TRUE. 
  LLPTENDR=.TRUE.
  LLPFRSO=.TRUE. 
  LLPFRTH=.TRUE. 
  LLPFPLSL=.TRUE. 
  LLPFPLSN=.TRUE. 
  LLPFPLSG=.TRUE.  
  LLPSTRTU=.TRUE. 
  LLPSTRTV=.TRUE. 
  LLPFRTHDS=.TRUE. 
  LLPFRSODS=.TRUE.  
  LLPCLCH=.TRUE. 
  LLPCLCL=.TRUE. 
  LLPCLCM=.TRUE. 
  LLPCLCT=.TRUE. 
  LLPFCS=.TRUE. 
  LLPFCLL=.TRUE. 
  LLPFCLN=.TRUE. 
  LLPFEVL=.TRUE. 
  LLPFEVN=.TRUE.
  
  ILEV(:)=0
  ILON(:)=0
  CLVAR(:)=' '
  IERRN=0
  
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      IF(ISNAN(PTENDT(JLON,JLEV)).AND.LLPTENDT.AND.IERRN < JPERRN) THEN
        LLPTENDT=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PTENDT'
      ENDIF
      IF(ISNAN(PTENDTKE(JLON,JLEV)).AND.LLPTENDTKE.AND.IERRN < JPERRN) THEN
        LLPTENDTKE=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PTENDTKE'
      ENDIF
      IF(ISNAN(PTENDR(JLON,JLEV,1)).OR.ISNAN(PTENDR(JLON,JLEV,2)).OR.ISNAN(PTENDR(JLON,JLEV,3)).OR.ISNAN(PTENDR(JLON,JLEV,4)).OR.&
        &ISNAN(PTENDR(JLON,JLEV,5)).OR.ISNAN(PTENDR(JLON,JLEV,6)).AND.LLPTENDR.AND.IERRN < JPERRN) THEN
        LLPTENDR=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PTENDR'
      ENDIF
    ENDDO
  ENDDO
  DO JLEV=0,KLEV
    DO JLON=KIDIA,KFDIA
      IF(ISNAN(PFRSO(JLON,JLEV)).AND.LLPFRSO.AND.IERRN < JPERRN) THEN
        LLPFRSO=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PFRSO'
      ENDIF
      IF(ISNAN(PFRTH(JLON,JLEV)).AND.LLPFRTH.AND.IERRN < JPERRN) THEN
        LLPFRTH=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PFRTH'
      ENDIF
      IF(ISNAN(PFPLSL(JLON,JLEV)).AND.LLPFPLSL.AND.IERRN < JPERRN) THEN
        LLPFPLSL=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PFPLSL'
      ENDIF
      IF(ISNAN(PFPLSN(JLON,JLEV)).AND.LLPFPLSN.AND.IERRN < JPERRN) THEN
        LLPFPLSN=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PFPLSN'
      ENDIF
      IF(ISNAN(PFPLSG(JLON,JLEV)).AND.LLPFPLSG.AND.IERRN < JPERRN) THEN
        LLPFPLSG=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PFPLSG'
      ENDIF
      IF(ISNAN(PSTRTU(JLON,JLEV)).AND.LLPSTRTU.AND.IERRN < JPERRN) THEN
        LLPSTRTU=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PSTRTU'
      ENDIF
      IF(ISNAN(PSTRTV(JLON,JLEV)).AND.LLPSTRTV.AND.IERRN < JPERRN) THEN
        LLPSTRTV=.FALSE.
        IERRN=IERRN+1
        ILEV(IERRN)=JLEV
        ILON(IERRN)=JLON
        CLVAR(IERRN)='PSTRTV'
      ENDIF
    ENDDO
  ENDDO

  DO JLON=KIDIA,KFDIA
    IF(ISNAN(PFRTHDS(JLON)).AND.LLPFRTHDS.AND.IERRN < JPERRN) THEN
      LLPFRTHDS=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PFRTHDS'
    ENDIF
    IF(ISNAN(PFRSODS(JLON)).AND.LLPFRSODS.AND.IERRN < JPERRN) THEN
      LLPFRSODS=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PFRSODS'
    ENDIF
    IF(ISNAN(PCLCH(JLON)).AND.LLPCLCH.AND.IERRN < JPERRN) THEN
      LLPCLCH=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PCLCH'
    ENDIF
    IF(ISNAN(PCLCL(JLON)).AND.LLPCLCL.AND.IERRN < JPERRN) THEN
      LLPCLCL=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PCLCL'
    ENDIF
    IF(ISNAN(PCLCM(JLON)).AND.LLPCLCM.AND.IERRN < JPERRN) THEN
      LLPCLCM=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PCLCM'
    ENDIF
    IF(ISNAN(PCLCT(JLON)).AND.LLPCLCT.AND.IERRN < JPERRN) THEN
      LLPCLCT=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PCLCT'
    ENDIF
    IF(ISNAN(PFCS(JLON)).AND.LLPFCS.AND.IERRN < JPERRN) THEN
      LLPFCS=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=JLEV
      ILON(IERRN)=KLON
      CLVAR(IERRN)='PFCS'
    ENDIF
    IF(ISNAN(PFCLL(JLON)).AND.LLPFCLL.AND.IERRN < JPERRN) THEN
      LLPFCLL=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PFCLL'
    ENDIF
    IF(ISNAN(PFCLN(JLON)).AND.LLPFCLN.AND.IERRN < JPERRN) THEN
      LLPFCLN=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PFCLN'
    ENDIF
    IF(ISNAN(PFEVL(JLON)).AND.LLPFEVL.AND.IERRN < JPERRN) THEN
      LLPFEVL=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PFEVL'
    ENDIF
    IF(ISNAN(PFEVN(JLON)).AND.LLPFEVN.AND.IERRN < JPERRN) THEN
      LLPFEVN=.FALSE.
      IERRN=IERRN+1
      ILEV(IERRN)=KLEV
      ILON(IERRN)=JLON
      CLVAR(IERRN)='PFEVN'
    ENDIF
  ENDDO
  
  DO JERRN=1,IERRN
    !
    !-------------------------------------------------
    ! A NaN has been found.
    !-------------------------------------------------
    !
    WRITE(NULERR,*) ' '
    WRITE(NULERR,*) 'CHECKNAN/ERROR : ',TRIM(CLVAR(JERRN)),' contains NaN value at ILEV=',ILEV(JERRN),' ILON=',ILON(JERRN)
    WRITE(NULERR,*) ' '
    WRITE(NULERR,*) 'Time of error: '
    WRITE(NULERR,*) '  NINDAT=',NINDAT
    WRITE(NULERR,*) '  KSTEP=',KSTEP
    WRITE(NULERR,*) '  TSPHY (physical time step)=',YDPHY2%TSPHY
    WRITE(NULERR,*) '  RSTATI (time in seconds since model integration start)=',RSTATI
    WRITE(NULERR,*) '  RHGMT (GMT time in seconds between 0 and 86400)=',RHGMT
    WRITE(NULERR,*) '  GMT time (in hours)=',RHGMT/3600._JPRB
    IF (PRESENT(ERRLOC)) THEN
      WRITE(NULERR,*) 'Error occurred at location: ', ERRLOC
    ENDIF
    !
    !-------------------------------------------------
    ! ZLOCST : local solar time in hours.
    ! The constant 0.261799387 rad is 15 degrees (the Earth rotates 15 degrees per hour).
    !-------------------------------------------------
    !
    ZLOCST=MODULO(RHGMT,86400._JPRB)/3600.+PGELAM(ILON(JERRN))/0.261799387_JPRB
    WRITE(NULERR,*) '  Local solar time (in hours)=',ZLOCST
    !
    !-------------------------------------------------
    ! Location of error.
    !-------------------------------------------------
    !
    WRITE(NULERR,*) ' '
    WRITE(NULERR,*) 'Vertical location of error: '
    WRITE(NULERR,*) '  ILEV=',ILEV(JERRN)
    WRITE(NULERR,*) '  KLEV=',KLEV
    WRITE(NULERR,*) '  PAPHIF(',ILON(JERRN),',',ILEV(JERRN),')=',PAPHIF(ILON(JERRN),ILEV(JERRN))
    WRITE(NULERR,*) '  Geopotential at surface : PAPHI(',ILON(JERRN),',KLEV)=',PAPHIF(ILON(JERRN),KLEV)
    WRITE(NULERR,*) '  Land-sea mask PLSM=',PLSM(ILON(JERRN))
    WRITE(NULERR,*) ' '
    WRITE(NULERR,*) 'Horizontal location of error: '
    WRITE(NULERR,*) '  KIDIA=',KIDIA
    WRITE(NULERR,*) '  KFDIA=',KFDIA
    WRITE(NULERR,*) '  ILON(JERRN)=',ILON(JERRN)
    ZCONRD=45._JPRB/ATAN(1._JPRB)
    ZLON=PGELAM(ILON(JERRN))*ZCONRD
    IF(ZLON > 180._JPRB) ZLON=ZLON-360._JPRB
    ZLAT=ZCONRD*ASIN(PGEMU(ILON(JERRN)))
    WRITE(NULERR,fmt='(2(A,F11.5))') '  Latitude, Longitude (in degrees) : ',ZLAT,' , ',ZLON
    WRITE(NULERR,*) '  PMU0 (cosine of Sun zenithal angle)=',PMU0(ILON(JERRN))
    !
    !-------------------------------------------------
    ! Vertical profile.
    !-------------------------------------------------
    !
    WRITE(NULERR,*) ' '
    WRITE(NULERR,*) 'Vertical profile of p, T and relative humidity : '
    DO JLEV=1,KLEV
      ZHR=PQ(ILON(JERRN),JLEV)/MAX(1.E-6_JPRB,PQSAT(ILON(JERRN),JLEV))
      WRITE(NULERR,FMT=*) PAPRSF(ILON(JERRN),JLEV),PT(ILON(JERRN),JLEV),ZHR
    ENDDO
    !
    !-------------------------------------------------
    ! Surface temperature.
    !-------------------------------------------------
    !
    WRITE(NULERR,*) ' '
    WRITE(NULERR,*) 'Surface temperature : ',PTS(ILON(JERRN))
  ENDDO
  !
  !-------------------------------------------------
  ! Generates abort.
  !-------------------------------------------------
  !
  CALL ABOR1('CHECKNAN : a NaN found in physical fluxes !')
ENDIF ! LLOK

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHECKNANARO',1,ZHOOK_HANDLE)
END SUBROUTINE CHECKNANARO
