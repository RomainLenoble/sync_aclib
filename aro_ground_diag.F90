SUBROUTINE ARO_GROUND_DIAG( KBL, KLON, KST, KEND, &
     & KLEV, KORDER, KDGUNG, KDGUXG,KDLUNG, KDLUXG, &
     & PZS, PFEV, PUM, PVM, PHEIGHT, PLWNET, PSWNET,PINDX, PINDY, &
     & PQS, PGZ0, PGZ0H, PT2M, PQ2M, PRH2M, PU10M, PV10M, PU10N, PV10N, &
     & PFCLL, PFCLN, PFEVL, PFEVN, PSSO_STDEV, PTWSNOW , PTDSNOW, PBUDLW, PBUDSW, &
     & PLE, PTOWN, PCD, PCH, PSIC )

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!!****  *AC_GROUND_PARAM*

!!    PURPOSE
!!    -------
!      Get surface diagnostics for ALADIN/ARPEGE

!!**  METHOD
!!    ------

!!    EXTERNAL
!!    --------

!!    IMPLICIT ARGUMENTS
!!    ------------------

!!    REFERENCE
!!    ---------

!!    AUTHOR
!!    ------
!!      L. Kraljevic          * Meteo-France *

!!    MODIFICATIONS
!!    -------------
!!      Original    08/2006
!       F. Bouyssel 05/2007 Use of PZS, PU, PV, PRH2M, PFCLL, PFCLN
!       M. Jidane   08/2008 Z0 and Z0H recovery from the various tiles
!       S. Riette   01/2010 Wind interpolation when 10m level is higher than
!                           forcing height
!       Y. Seity    03/2010 add PFEVN and PFEVL
!       F. Bouyssel 05/2010 Use of PINDX, PINDY
!       P.Marguinaud 10/2010 OpenMP + cleaning
!       S. Riette   06/2010 PSSO_STDEV and PTWSNOW added
!       B. Decharme 12/2009 Atmospheric orography interpolation of T2M switched
!       B. Decharme 12/2009 Correction of the atmospheric net thermic flux
!       B. Decharme 12/2010 Pass latent heat flux to atmos (not recalculated in aplpar)
!       I. Etchevers 10/2022 PTDSNOW, snow depth, added
!       J. Masek    08/2023 Removed linear averaging of roughness lengths over
!                           tiles; averages from GET_SURF_VAR_N are taken
!-------------------------------------------------------------------------------

!*       0.     DECLARATIONS
!               ------------

USE MODD_SURF_PAR,     ONLY :  XUNDEF
USE MODD_CSTS,         ONLY : XG, XTT, XKARMAN
USE MODD_PREP,         ONLY : XT_CLIM_GRAD
USE MODD_FROMMPA,      ONLY : CPROGRAM_MSE
USE MODD_SURFEX_ARO,   ONLY : YSURFEX_ARO_ALL, YSURFEX_ARO_CUR

USE MODI_GET_SURF_VAR_N
USE MODI_GET_FLUX_N
USE MODI_DIAG_SURF_ATM_N

USE YOMCT0, ONLY : LELAM

IMPLICIT NONE

!*      0.1    declarations of arguments


INTEGER(KIND=JPIM),                 INTENT(IN)  :: KBL       ! NPROMA-packat number
INTEGER(KIND=JPIM),                 INTENT(IN)  :: KLON      ! horizontal dimension
INTEGER(KIND=JPIM),                 INTENT(IN)  :: KST       ! start of physical values
INTEGER(KIND=JPIM),                 INTENT(IN)  :: KEND      ! end of physical values
INTEGER(KIND=JPIM),                 INTENT(IN)  :: KLEV      ! number of vertcal levels given
INTEGER(KIND=JPIM),                 INTENT(IN)  :: KORDER    ! 1 when levels are from
                                                  ! bottom to top, -1 otherwise
INTEGER(KIND=JPIM),                 INTENT(IN)  :: KDGUNG, KDGUXG,KDLUNG, KDLUXG !dim param
                                                  ! to fill E-Zone
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(IN)  :: PZS       ! orography (of atmospheric model)
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(IN)  :: PFEV      ! surface flux of water vapor (m/s*kg/kg)
REAL(KIND=JPRB), DIMENSION(KLON,KLEV),   INTENT(IN)  :: PUM  ! zonal wind at lowest model level
REAL(KIND=JPRB), DIMENSION(KLON,KLEV),   INTENT(IN)  :: PVM  ! meridional wind at lowest model level
REAL(KIND=JPRB), DIMENSION(KLON,KLEV),   INTENT(IN)  :: PHEIGHT ! height of U and V levels
REAL(KIND=JPRB), DIMENSION(KLON),        INTENT(IN)  :: PLWNET  ! atmospheric net thermic flux
REAL(KIND=JPRB), DIMENSION(KLON),        INTENT(IN)  :: PSWNET  ! atmospheric net shortwave flux
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(IN)  :: PINDX     ! x coordinate on the plane
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(IN)  :: PINDY     ! y coordinate on the plane
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PQS       ! surface humidity
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PGZ0      ! rougness lenght for momentum
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PGZ0H     ! rouhness length for heat
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PT2M      ! temperature at 2m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PQ2M      ! specific humidity at 2m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PRH2M     ! relative humidity at 2m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PU10M     ! zonal wind at 10m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PV10M     ! meridional wind at 10m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PU10N     ! zonal neutral wind at 10m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PV10N     ! meridional neutral wind at 10m
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PFCLL     ! flux de chaleur latente (liquide)
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PFCLN     ! flux de chaleur latente (solide)
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PFEVL     ! flux de vapeur d'eau sur eau liq.
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PFEVN     ! flux de vapeur d'eau sur neige
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PSSO_STDEV! S.S.O. standard deviation (m)
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PTWSNOW   ! Snow total reservoir
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PTDSNOW   ! Snow depth
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PSIC      ! Sea ice concentration
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PBUDLW    ! Longwave budget betwwen ARP and SURFEX
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PBUDSW    ! Shortwave budget betwwen ARP and SURFEX
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PLE       ! Surface latent heat flux
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PTOWN     ! Fraction that is town
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PCD       ! drag coeficient at the surface
REAL(KIND=JPRB), DIMENSION(KLON),   INTENT(OUT) :: PCH       ! heat coeficient at the surface

!-------------------------------------------------------------------------------

!*      0.2    declarations of local variables
!              -------------------------------
INTEGER :: JI   ! loop counter
INTEGER :: ILEVEL ! loop counter for levels
INTEGER :: IBOT, ITOP !bottom and top model levels

!* Atmospheric variables
!  ---------------------

!            suffix 'A' stands for atmospheric variable at first model level

!* Arrays put in 1D vectors
!  ------------------------

REAL, DIMENSION(KLON)    :: ZZS
REAL, DIMENSION(KLON)    :: ZQS_SEA
REAL, DIMENSION(KLON)    :: ZQS_WATER
REAL, DIMENSION(KLON)    :: ZQS_NATURE
REAL, DIMENSION(KLON)    :: ZQS_TOWN
REAL, DIMENSION(KLON)    :: ZSEA
REAL, DIMENSION(KLON)    :: ZWATER
REAL, DIMENSION(KLON)    :: ZNATURE
REAL, DIMENSION(KLON)    :: ZTOWN
REAL, DIMENSION(KLON)    :: ZRN
REAL, DIMENSION(KLON)    :: ZH
REAL, DIMENSION(KLON)    :: ZLE
REAL, DIMENSION(KLON)    :: ZLEI
REAL, DIMENSION(KLON)    :: ZGFLUX
REAL, DIMENSION(KLON)    :: ZDELT
REAL, DIMENSION(KLON)    :: ZU10M
REAL, DIMENSION(KLON)    :: ZV10M
REAL, DIMENSION(KLON)    :: ZVMODM
REAL, DIMENSION(KLON)    :: ZVMOD10
REAL, DIMENSION(KLON)    :: ZSURFLWNET
REAL, DIMENSION(KLON)    :: ZSURFSWNET
REAL, DIMENSION(KLON)    :: ZEVAP
REAL, DIMENSION(KLON)    :: ZSUBL
REAL, DIMENSION(KLON)    :: ZIV



!* Other real variables
!  -----------------------
REAL :: ZDGUN            ! y dimension of C+I+E zone
REAL :: ZDGUX            ! y dimension of C+I zone
REAL :: ZDLUN            ! x dimension of C+I+E zone
REAL :: ZDLUX            ! x dimension of C+I zone

REAL, DIMENSION(KLON)    :: ZEMASK  ! E-zone mask

!-------------------------------------------------------------------------------

!*CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_GROUND_DIAG',0,ZHOOK_HANDLE)

PBUDLW(:) = 0.0
PBUDSW(:) = 0.0

YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (KBL)

CALL DIAG_SURF_ATM_N &
  (YSURFEX_ARO_CUR, &
   CPROGRAM_MSE)
!  (YSURFEX_ARO_CUR%IM%DGEI, &
!   YSURFEX_ARO_CUR%FM%DGF,  &
!   YSURFEX_ARO_CUR%DGL,     &
!   YSURFEX_ARO_CUR%IM%DGI,  &
!   YSURFEX_ARO_CUR%SM%DGS,  &
!   YSURFEX_ARO_CUR%DGU,     &
!   YSURFEX_ARO_CUR%TM%DGT,  &
!   YSURFEX_ARO_CUR%WM%DGW,  &
!   YSURFEX_ARO_CUR%U,       &
!   YSURFEX_ARO_CUR%USS,     &
!   CPROGRAM_MSE)

CALL GET_SURF_VAR_N &
  (YSURFEX_ARO_CUR%FM, YSURFEX_ARO_CUR%IM,    &
   YSURFEX_ARO_CUR%SM, YSURFEX_ARO_CUR%TM, &
   YSURFEX_ARO_CUR%WM, YSURFEX_ARO_CUR%DUO,     &
   YSURFEX_ARO_CUR%DU, YSURFEX_ARO_CUR%UG,  &
   YSURFEX_ARO_CUR%U,      YSURFEX_ARO_CUR%USS,     &
   YSURFEX_ARO_CUR%SM%S,                            &
   CPROGRAM_MSE,KLON,1, PZS=ZZS,                    &
   PZ0=PGZ0,PZ0H=PGZ0H,                             &
   PQS_SEA    =ZQS_SEA,     PQS_WATER =ZQS_WATER,   &
   PQS_NATURE =ZQS_NATURE,  PQS_TOWN  =ZQS_TOWN,    &
   PSEA=ZSEA, PWATER=ZWATER,                        &
   PNATURE=ZNATURE, PTOWN=ZTOWN, PSSO_STDEV=PSSO_STDEV,&
   PTWSNOW=PTWSNOW, PTDSNOW=PTDSNOW, PSIC=PSIC)

CALL GET_FLUX_N &
   (YSURFEX_ARO_CUR%DUO, YSURFEX_ARO_CUR%DU, CPROGRAM_MSE, KLON, ZRN, ZH,         &
   ZLE, ZLEI, ZGFLUX, PT2M, PQ2M, PRH2M, ZU10M,               &
   ZV10M,PSURFLWNET=ZSURFLWNET,PSURFSWNET=ZSURFSWNET,PCD=PCD, &
   PCH=PCH,PEVAP=ZEVAP,PSUBL=ZSUBL)

PLE= -1.0 * ZLE

PTOWN = ZTOWN

! #############################################################################
! Au lieu de faire juste une moyenne lineaire de ces champs,
! il vaudrait mieux, pour QS, passer par une moyenne des flux de chaleur et faire
! la derive ensuite ...

PQS = ZSEA   *ZQS_SEA    + &
    & ZWATER *ZQS_WATER  + &
    & ZNATURE*ZQS_NATURE + &
    & ZTOWN  *ZQS_TOWN

    ! print*, "ARO_GROUND_DIAG", PQS, ZSEA, ZQS_SEA, ZWATER , ZQS_WATER, ZNATURE, ZQS_NATURE, ZTOWN, ZQS_TOWN

! #############################################################################

PT2M = PT2M + &
   (PZS - ZZS) * XT_CLIM_GRAD

!  When 10m is higher than the forcing level, SURFEX sends the XUNDEF value for
!  ZU10M and ZV10M. Where these values are found, starting from the bottom, we
!  look for the model level just below 10m and we interpolate linearly the
!  model wind field.

!Bottom and top levels
IF(KORDER==1) THEN
  IBOT=1
  ITOP=KLEV
ELSE !must be KORDER=-1
  IBOT=KLEV
  ITOP=1
ENDIF

ILEVEL=IBOT
!While there are XUNDEF values and we aren't at model's top
DO WHILE(ILEVEL/=ITOP)

  !Where interpolation is needed and possible
  !(10m is between ILEVEL and ILEVEL+KORDER or 10m is below the bottom level)
  WHERE(PHEIGHT(:,ILEVEL)<=10._JPRB .AND. &
                                PHEIGHT(:,ILEVEL+KORDER)>10._JPRB)

    !Interpolation between ILEVEL and ILEVEL+KORDER
    ZU10M=PUM(:,ILEVEL) + &
          (PUM(:,ILEVEL+KORDER)-PUM(:,ILEVEL)) * &
          (10._JPRB-PHEIGHT(:,ILEVEL)) / &
          (PHEIGHT(:,ILEVEL+KORDER)-PHEIGHT(:,ILEVEL))
    ZV10M=PVM(:,ILEVEL) + &
          (PVM(:,ILEVEL+KORDER)-PVM(:,ILEVEL)) * &
          (10._JPRB-PHEIGHT(:,ILEVEL)) / &
          (PHEIGHT(:,ILEVEL+KORDER)-PHEIGHT(:,ILEVEL))
  END WHERE
  ILEVEL=ILEVEL+KORDER !level just higher
END DO

ZVMOD10=SQRT(ZU10M**2 + ZV10M**2)
ZVMODM = &
          SQRT(PUM(:,IBOT)**2 + PVM(:,IBOT)**2)
WHERE (ZVMODM>0._JPRB)
  PU10M = PUM(:,IBOT) * &
                         ZVMOD10 / ZVMODM
  PV10M = PVM(:,IBOT) * &
                         ZVMOD10 / ZVMODM
ELSE WHERE
  PU10M = 0._JPRB
  PV10M = 0._JPRB
END WHERE


PGZ0  = XG * PGZ0
PGZ0H = XG * PGZ0H

! Radiative budget between ARP and SURFEX
PBUDLW = ZSURFLWNET - PLWNET
PBUDSW = ZSURFSWNET - PSWNET


ZDELT = MAX(0.0_JPRB, SIGN(1.0,XTT-PT2M) )
PFCLL = - ZLE * ( 1.0_JPRB - ZDELT )
PFCLN = - ZLE * ZDELT
PFEVL = - PFEV * ( 1.0_JPRB - ZDELT )
PFEVN = - PFEV * ZDELT

IF (LELAM) THEN
!  In E zone

  ZDGUN = REAL(KDGUNG)
  ZDGUX = REAL(KDGUXG)
  ZDLUN = REAL(KDLUNG)
  ZDLUX = REAL(KDLUXG)

  DO JI=KST,KEND
    ZEMASK(JI)=MAX(0.0_JPRB,SIGN(1.0,(PINDX(JI)-ZDLUN)*(ZDLUX-PINDX(JI))))&
     & *MAX(0.0_JPRB,SIGN(1.0,(PINDY(JI)-ZDGUN)*(ZDGUX-PINDY(JI))))
  ENDDO

  WHERE (ZEMASK(:) == 0._JPRB)
    PQS  (:) = 0.0_JPRB
    PGZ0 (:) = 1.0_JPRB
    PGZ0H(:) = 0.1_JPRB
    PT2M (:) = XTT
    PQ2M (:) = 0.0_JPRB
    PRH2M(:) = 0.0_JPRB
    PU10M(:) = 0.0_JPRB
    PV10M(:) = 0.0_JPRB
    PFCLL(:) = 0.0_JPRB
    PFCLN(:) = 0.0_JPRB
    PSSO_STDEV(:) = 0.0_JPRB
    PTWSNOW(:) = 0.0_JPRB
    PTDSNOW(:) = 0.0_JPRB
    PSIC(:) = 0.0_JPRB
    PBUDLW(:) = 0.0_JPRB
    PBUDSW(:) = 0.0_JPRB
    PLE(:) = 0.0_JPRB
    PCD(:) = 0.0_JPRB
  END WHERE
ENDIF

WHERE (PTWSNOW(:) == XUNDEF)
  PTWSNOW(:) = 0.0_JPRB
ENDWHERE

WHERE (PTDSNOW(:) == XUNDEF)
  PTDSNOW(:) = 0.0_JPRB
ENDWHERE

WHERE (PSIC(:) == XUNDEF)
  PSIC(:) = 0.0_JPRB
ENDWHERE

! Neutral wind computation
!CP For scatterometry purposes (cy42_op)
!CP Must be moved on SURFEX routines later (PARAM_CLS/CLS_WIND routines), cy43?
!CP Computation here with Z0, instead of ZBN
!CP ZIV result may be higher than 1 (model bottom height lower than 10m)
!CP test ZIV(:)=MAX(0._JPRB,MIN(1.0_JPRB,LOG(1.+10._JPRB*XG/PGZ0(:))*SQRT(PCD(:))/XKARMAN))
!CP while this computation will be done here, that would be done after the E
! zone reset in E zone (LELAM case), else PGZ0 "undefined"
ZIV(:)=MAX(0._JPRB,LOG(1._JPRB+10._JPRB*XG/PGZ0(:))*SQRT(PCD(:))/XKARMAN)

!Same trick as for real wind
!CP warning: ZVMOD10 reused for the neutral case
ZVMOD10=SQRT((PUM(:,IBOT)*ZIV(:))**2 + (PVM(:,IBOT)*ZIV(:))**2)
!CP ZVMODM: already done for real wind

WHERE (ZVMODM>0._JPRB)
  PU10N = PUM(:,IBOT) * &
                         ZVMOD10 / ZVMODM
  PV10N = PVM(:,IBOT) * &
                         ZVMOD10 / ZVMODM
ELSE WHERE
  PU10N = 0._JPRB
  PV10N = 0._JPRB
END WHERE

IF (LELAM) THEN
!  In E zone

  !ZEMASK already defined
  WHERE (ZEMASK(:) == 0._JPRB)
    PU10N(:) = 0.0_JPRB
    PV10N(:) = 0.0_JPRB
  END WHERE
ENDIF
! End neutral wind computation

! Rotate the wind components (not done now)


IF (LHOOK) CALL DR_HOOK('ARO_GROUND_DIAG',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_GROUND_DIAG
