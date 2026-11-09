SUBROUTINE SUPHMSE(YDGEOMETRY,YDMODEL,KULOUT)

!**** *SUPHMSE*   - Initialize common meso_NH MODD_ used in physics for AROME

!     Purpose.
!     --------
!           Initialize MODD_PARAMETERS, MODD_CST, MODD_CONF,
!           MODD_RAIN_ICE_DESCR, MODD_RAIN_ICE_PARAM, MODD_BUDGET  
!           parameters used in meso_NH Physics and aladin/meso_NH physics
!           interface 

!**   Interface.
!     ----------
!        *CALL* *SUPHMSE(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!        R. Zaaboul
!        Original : 28-Feb-2006

!     Modifications.
!     --------------
!        Y. Seity : 28-March-2007 Move chemistry setup from suphmpa
!        A. Alias : 13-June-2007 Set IIEXT and IJEXT for LELAM=.F. and 
!                                new call to aroini_surf when ARPEGE
!        A. Alias : 26-Sept-2007 Call to aroini_surf modified
!        Y. Seity : 11-Jan-2008 Call to aroini_surf modified to check date
!        Y. Seity : 30-01-2008 Add arguments to aroini_surf 
!        M. Mkhtari:01-02-2011 add aroini_wet_dep and the key LMDUST   
!        B. Decharme : 12-2010 Add arguments to aroini_surf for Earth System Model 
!        B. Decharme : 12-2009 Zenithal angle as in ARPEGE/ALADIN 
!        R. El Khatib 16-Jun-2010 Namelist
!        P.Marguinaud 10-Aug-2010 More namelist parameters
!        Y. Seity : 14-Feb-2011 Init modd_frommpa and new args to aroini_surf
!        K. Yessad: Sep 2010 : organigramme simplification
!        A. Alias : Mar-2011 LASTRF to prevent any drift in insolation (A.Voldoire)
!        M. Jerczynski : Jun 2011 some cleaning to meet norms
!        P. Marguinaud : Jul-2011 Fix KSURFEXCTL option and add more parameters to NAMPHMSE
!        K. Yessad (July 2014): Move some variables.
!        2016-09, M. Mokhtari & A. Ambar: call for the routine
!                                         aroini_wet_dep.F90
!      R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!        2023-06, M. Michou  : consider prog aerosols and/or prog gazeous chemical species as tracers
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    ,ONLY : YDCST=>YRCST ! allows use of included functions. REK.
USE YOMRIP0   ,ONLY : NINDAT
USE YOMCT0    ,ONLY : LELAM, L_OOPS
USE YOMNSV    ,ONLY : NSV_CHEMBEG, NSV_CHEMEND, NSV_DSTBEG, NSV_DSTEND,&
 &                    NSV_AERBEG, NSV_AEREND, NSV_CO2,NSV_DSTDEPBEG,&
 &                    NSV_DSTDEPEND
USE YOMMP0    ,ONLY : MYPROC
USE MODE_INI_CST, ONLY: INI_CST
USE MODD_NSV, ONLY: NSV_ASSOCIATE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE (MODEL)      ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDAY, IMONTH, IYEAR, IMYPROC
INTEGER(KIND=JPIM) :: INSV_CHEMBEG,INSV_CHEMEND, INSV

LOGICAL :: LLNOTMAP
LOGICAL :: LLCHEMAERO_ARPC, LLCHEM_ARPC, LLAERO_ARPC 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "fctast.func.h"
#include "fcttim.func.h"
#include "aroini_mnhc.h"
#include "aroini_nsv.h"
#include "suphmse_surface.h"
#include "aroini_wet_dep.h"
#include "aroini_frommpa.h"
#include "aroini_conf.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPHMSE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEGEO=>YDGEOMETRY%YREGEO, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)
ASSOCIATE(LRDUST=>YDARPHY%LRDUST, LINITDUST=>YDARPHY%LINITDUST, &
 & LRDEPOS=>YDARPHY%LRDEPOS, LUSECHEM=>YDARPHY%LUSECHEM, &
 & LINITORILAM=>YDARPHY%LINITORILAM, LMDUST=>YDARPHY%LMDUST, &
 & LINITCHEM=>YDARPHY%LINITCHEM, LORILAM=>YDARPHY%LORILAM, LMPA=>YDARPHY%LMPA, &
 & NPROMA=>YDDIM%NPROMA, &
 & NGFL_EXT=>YGFL%NGFL_EXT, &
 & RDAY=>YDCST%RDAY, REA=>YDCST%REA, REPSM=>YDCST%REPSM, &
 & LRAYFM15=>YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM15, NLIMA=>YGFL%NLIMA, & 
 & LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM,NCHEM=>YDMODEL%YRML_GCONF%YGFL%NCHEM,CHEM_SCHEME=>YDMODEL%YRML_CHEM%YRCHEM%CHEM_SCHEME, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME,NAERO=>YGFL%NAERO)
!     ------------------------------------------------------------------

IF (LRAYFM15) CALL ABOR1('SUPHMSE: FM15 NOT SUPPORTED WITH SURFEX!')

! Initialize MNH constants if not LMPA

IF (LELAM) THEN
  LLNOTMAP=.NOT.YDGEOMETRY%YREGEO%LMAP
ELSE
  LLNOTMAP=.TRUE.
ENDIF
IF (.NOT.LMPA) THEN
  CALL INI_CST
ENDIF
CALL AROINI_CONF(KULOUT,YDDYNA%LTWOTL,LLNOTMAP)

!     initialisation of chemistry, aerosols and dust scheme
! ANNEE

WRITE(KULOUT,*)'NINDAT =',NINDAT
IYEAR = NINDAT / 10000
! MOIS

IMONTH = (NINDAT - 10000*IYEAR ) / 100
! JOUR DU MOIS

IDAY = NINDAT - 10000*IYEAR - 100*IMONTH

 IMYPROC=MYPROC

 CALL AROINI_MNHC(LUSECHEM, LORILAM, LRDUST, LRDEPOS,&
  &                LINITCHEM, LINITDUST, LINITORILAM,&
  &                IDAY, IMONTH, IYEAR, KULOUT,IMYPROC)

!       Initialisation of aerosols dust wet deposition for Aladin

IF (LMDUST.AND.(NGFL_EXT/=0).AND.LRDEPOS) THEN
  CALL AROINI_WET_DEP 
ENDIF

!       initialisation nsv
!       initialisation for when ARPEGE-Climat prognostic chemistry and/or prognostic aerosol are active
LLAERO_ARPC     = .FALSE.
LLCHEM_ARPC     = .FALSE.
LLCHEMAERO_ARPC = .FALSE.
INSV            = 0_JPIM

LLAERO_ARPC     = AERO_SCHEME == "tactic" .AND. NAERO>0
LLCHEM_ARPC     = LCHEM_ARPCLIM.AND.(TRIM(CHEM_SCHEME).EQ.'arpclim_relacs'.OR.TRIM(CHEM_SCHEME).EQ.'mocage')
LLCHEMAERO_ARPC = LLAERO_ARPC.AND.LLCHEM_ARPC 

IF (LLCHEMAERO_ARPC.OR.LLCHEM_ARPC.OR.LLAERO_ARPC) THEN ! case with interactive ARPEGE-Climat prognostic chemistry/aerosols and DUST scheme/other chemistry not coded ; case with CO2 (NSV_CO2) to be seen

    IF (LLCHEM_ARPC) THEN
       INSV_CHEMBEG = 1_JPIM
       INSV_CHEMEND = NCHEM
    ENDIF 
    IF (LLAERO_ARPC) THEN
       INSV_CHEMBEG = 1_JPIM
       INSV_CHEMEND = NAERO
    ENDIF    
    IF (LLCHEMAERO_ARPC) THEN
       INSV_CHEMBEG = 1_JPIM
       INSV_CHEMEND = NCHEM + NAERO
    ENDIF     

    CALL NSV_ASSOCIATE()
    CALL AROINI_NSV(NLIMA,INSV_CHEMBEG,INSV_CHEMEND,INSV,INSV, &
 &                   INSV,INSV,INSV,INSV,INSV,YDMODEL)
    WRITE(UNIT=KULOUT,FMT='('' NLIMA = '',I3,'' NSV_CHEMBEG = '',I3,&
 & '' NSV_CHEMEND = '',I3,'' NSV_AERBEG = '',I3,'' NSV_AEREND = '',I3,&
 & '' NSV_DSTBEG = '',I3,'' NSV_DSTEND = '',I3,&
 & '' NSV_DSTDEPBEG = '',I3,'' NSV_DSTDEPEND = '',I3,'' NSV_CO2 = '',I3)')&
 & NLIMA,NSV_CHEMBEG,NSV_CHEMEND,NSV_AERBEG,NSV_AEREND,NSV_DSTBEG,NSV_DSTEND,&
 & NSV_DSTDEPBEG,NSV_DSTDEPEND,NSV_CO2

ELSE 
     CALL NSV_ASSOCIATE()
     CALL AROINI_NSV(NLIMA,NSV_CHEMBEG, NSV_CHEMEND, NSV_AERBEG, NSV_AEREND,&
   &             NSV_DSTBEG, NSV_DSTEND, NSV_DSTDEPBEG, NSV_DSTDEPEND,&
   &             NSV_CO2,YDMODEL)

     WRITE(UNIT=KULOUT,FMT='('' NLIMA = '',I3,'' NSV_CHEMBEG = '',I3,&
 & '' NSV_CHEMEND = '',I3,'' NSV_AERBEG = '',I3,'' NSV_AEREND = '',I3,&
 & '' NSV_DSTBEG = '',I3,'' NSV_DSTEND = '',I3,&
 & '' NSV_DSTDEPBEG = '',I3,'' NSV_DSTDEPEND = '',I3,'' NSV_CO2 = '',I3)')&
 & NLIMA,NSV_CHEMBEG,NSV_CHEMEND,NSV_AERBEG,NSV_AEREND,NSV_DSTBEG,NSV_DSTEND,&
 & NSV_DSTDEPBEG,NSV_DSTDEPEND,NSV_CO2

ENDIF

!     Initialisation of variables from modd_frommpa.F90
CALL AROINI_FROMMPA

!     Surface setup ( Moved into MODEL_INIT in the OOPS case )
IF (.NOT.L_OOPS) THEN
  CALL SUPHMSE_SURFACE(YDGEOMETRY,YDMODEL,KULOUT, 'C', NPROMA)
ENDIF


! -----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHMSE',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHMSE
