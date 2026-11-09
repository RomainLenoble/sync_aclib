SUBROUTINE AER_MAIN1(YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS, &
            & KIDIA, KFDIA, &
            & KLUOUT, &
            & KLON, KLEV, PTSTEP, &
            & PAP, PAPH, PAPHI, PDELP, PTP, PVERVEL, PRH, PRHODREFM, &
            & PSST, PWS1, PCI, PBARESOIL, PLSM, PCLK, PSNS, &
            & PAERDEP, PAERSCC, PZ0M, PWIND, &
            & PAER, PTENAER, PSFAER, &
            & YAEROUT_ACLIB )

! INPUTS:
! -------
! KIDIA , KFDIA           : BEGIN/END of LON BLOCK
! KLUOUT                  : UNIT OF OUTPUT FILE
! KLON                    : NUMBER OF LONGITUDES / LATITUDES
! KLEV                    : NUMBER OF VERTICAL LEVELS
! PTSTEP                  : TIME STEP 
! PAP (KLON,KLEV)         : FULL-LEVEL PRESSURE           (Pa)
! PAPH (KLON,0:KLEV)      : HALF-LEVEL PRESSURE           (Pa)
! PAPHI(KLON,0:KLEV)      : GEOPOTENTIAL ON HALF-LEVELS 
! PDELP (KLON,KLEV)       : DELTA PRESSURE                (Pa)
! PTP(KLON,KLEV)          : FULL-LEVEL TEMPERATURE        (K)
! PVERVEL(KLON,KLEV)      : FULL-LEVEL VERTICAL VELOCITY  (Pa s-1)
! PRH(KLON,KLEV)          : RELATIVE HUMIDITY
! PRHODREFM(KLON,KLEV)    : Dry air density
! PSST                    : Sea surface temperature (K)
! PWS1                    : Soil humidity
! PCI(KLON)               : FRACTION OF SEA-ICE
! PBARESOIL(KLON)         : BARE SOIL FRACTION
! PLSM(KLON)              : LAND-SEA MASK (note: treats lakes as sea)
! PCLK(KLON)              : LAKE COVER
! PSNS(KLON)              : SNOW DEPTH
! PAERDEP(KLON)           : SAND CONTENT (%)
! PAERSCC(KLON)           : CLAY CONTENT (%)
! PZ0M(KLON)              : ROUGHNESS LENGTH (m)
! PWIND(KLON)             : SURFACE WIND SPEED           (m/s)   

! INPUT/OUTPUT:
! -------
! PAER(KLON,KLEV,NAERO)    : AEROSOL CONCENTRATION
! PTENAER(KLON,KLEV,NAERO) : AEROSOL TENDENCY
! PSFAER(KLON,NAERO)       : AEROSOL SURFACE FLUX
! YAEROUT_ACLIB            : STRUCTURE FOR AEROSOL OUTPUT


USE PARKIND1,             ONLY : JPIM     ,JPRB
USE YOMHOOK,              ONLY : LHOOK,   DR_HOOK
USE M_AER_ACLIB,          ONLY : AER_ACLIB
USE MODD_CST_ACLIB,       ONLY : TCST_ACLIB
USE MODD_PARAM_AER_ACLIB, ONLY : TPARAM_AER_ACLIB, TAEROSOL, TAEROUT_ACLIB
USE, INTRINSIC :: IEEE_ARITHMETIC

IMPLICIT NONE

TYPE(TCST_ACLIB),        INTENT(IN)  :: YDCST_ACLIB
TYPE(TPARAM_AER_ACLIB),  INTENT(IN)  :: YPAER_ACLIB
TYPE(TAEROSOL),DIMENSION(YPAER_ACLIB%NAERO), INTENT(IN)  :: YAEROSOLS
INTEGER(KIND=JPIM),      INTENT(IN)  :: KIDIA , KFDIA 
INTEGER(KIND=JPIM),      INTENT(IN)  :: KLUOUT
INTEGER(KIND=JPIM),      INTENT(IN)  :: KLON , KLEV
REAL(KIND=JPRB),         INTENT(IN)  :: PTSTEP
REAL(KIND=JPRB),         INTENT(IN)  :: PAP(KFDIA,KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PAPH(KFDIA,0:KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PAPHI(KFDIA,0:KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PDELP(KFDIA,KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PTP(KFDIA,KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PVERVEL(KFDIA,KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PRH(KFDIA,KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PRHODREFM(KFDIA,KLEV)   
REAL(KIND=JPRB),         INTENT(IN)  :: PSST(KFDIA) 
REAL(KIND=JPRB),         INTENT(IN)  :: PWS1(KFDIA) 
REAL(KIND=JPRB),         INTENT(IN)  :: PCI(KFDIA) 
REAL(KIND=JPRB),         INTENT(IN)  :: PBARESOIL(KFDIA) 
REAL(KIND=JPRB),         INTENT(IN)  :: PLSM(KFDIA) 
REAL(KIND=JPRB),         INTENT(IN)  :: PCLK(KFDIA) 
REAL(KIND=JPRB),         INTENT(IN)  :: PSNS(KFDIA)
REAL(KIND=JPRB),         INTENT(IN)  :: PAERDEP(KFDIA)
REAL(KIND=JPRB),         INTENT(IN)  :: PAERSCC(KFDIA)
REAL(KIND=JPRB),         INTENT(IN)  :: PZ0M(KFDIA)
REAL(KIND=JPRB),         INTENT(IN)  :: PWIND(KFDIA)

REAL(KIND=JPRB),         INTENT(INOUT)  :: PAER(KFDIA,KLEV,YPAER_ACLIB%NAERO)   
REAL(KIND=JPRB),         INTENT(INOUT)  :: PTENAER(KFDIA,KLEV,YPAER_ACLIB%NAERO)   
REAL(KIND=JPRB),         INTENT(INOUT)  :: PSFAER(KFDIA,YPAER_ACLIB%NAERO)   
TYPE(TAEROUT_ACLIB),     INTENT(INOUT)  :: YAEROUT_ACLIB


REAL(KIND=JPRB) :: ZHOOK_HANDLE, val

! Local variables
INTEGER(KIND=JPIM)           :: JLEV,JLON,JAER
INTEGER(KIND=JPIM)           :: IAER_S, IAER_E, INAER
INTEGER(KIND=JPIM)           :: IDUMMY
REAL(KIND=JPRB)              :: ZRHO(KLON,KLEV)
REAL(KIND=JPRB)              :: ZDUMMY1(KLON,KLEV)
 
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_MAIN1',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!
!! Initialisation !!
!!!!!!!!!!!!!!!!!!!!

  !  DO JLEV=1,KLEV
  !    DO JLON=KIDIA,KFDIA
  !    !   ZRHO(JLON,JLEV)=PAP(JLON,JLEV)/(YDCST_ACLIB%RRD*PTP(JLON,JLEV))
  !      ZRHO(JLON,JLEV)=PRHODREFM(JLON, JLEV)
  !         IF (ZRHO(JLON,JLEV) > HUGE(1._JPRB)) THEN
  !         print*, "AERMain1, ZRHO infinity", ZRHO(JLON,JLEV), PAP(JLON,JLEV), YDCST_ACLIB%RRD, PTP(JLON,JLEV), JLON, JLEV
  !      ENDIF
  !    ENDDO
  !  ENDDO

 DO JLEV=1,KLEV
   DO JLON=KIDIA,KFDIA
        IF (1/(PAPH(JLON,JLEV-1)-PAPH(JLON,JLEV)) > HUGE(1._JPRB)) THEN
        print*, "AERMain1, PAPH issue ", JLON, JLEV, PAPH(JLON,:) 
     ENDIF
   ENDDO
ENDDO
!!!!!!!!!!!!!!!!!!
!! Dust family  !!
!!!!!!!!!!!!!!!!!!


DO JLON=KIDIA,KFDIA
  DO JAER=1, YPAER_ACLIB%NAERO

      val = PSFAER(JLON,JAER)

      ! Check NaN or Inf
      IF (.NOT. IEEE_IS_FINITE(val)) THEN
        PRINT *, "AerMain1 NaN/Inf detected in PSFAER at before aclib:", JLON, JAER, val
      END IF

      ! Check absurd values (adjust threshold!)
      IF (ABS(val) > 1.E10) THEN
        PRINT *, "AerMain1 Huge value in PSFAER at before aclib:", JLON, JAER, val
      END IF

    END DO
  END DO


   IF (YPAER_ACLIB%LDUSTON) THEN
       IAER_S = YPAER_ACLIB%NDUST_F0
       IAER_E = YPAER_ACLIB%NDUST_F1
       INAER = IAER_E - IAER_S + 1

       ! Emissions
      !  WRITE(KLUOUT,*) "ACLIB : CALL DUST EMIS (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( & 
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "EMIS", YPAER_ACLIB%CDUST_EMIS_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, &
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PRH=PRH, PRHODREF=PRHODREFM, &
            & PWIND=PWIND, PSST=PSST, PWS1=PWS1, PCI=PCI, &
            & PBARESOIL=PBARESOIL, PLSM=PLSM, PCLK=PCLK, PSNS=PSNS, &
            & PAERDEP=PAERDEP, PAERSCC=PAERSCC, PZ0M=PZ0M, &
            & PSFAER=PSFAER(:,IAER_S:IAER_E) &
            & )

       ! Dry deposition
      !  WRITE(KLUOUT,*) "ACLIB : CALL DUST DRYDEP (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( & 
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "DRYDEP", YPAER_ACLIB%CDUST_DRYDEP_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, &
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PAP=PAP, PAPH=PAPH, PAPHI=PAPHI, PDELP=PDELP, &
            & PRHODREF=PRHODREFM, &
            & PCI=PCI, PLSM=PLSM, PSNS=PSNS, &
            & PSFAER=PSFAER(:,IAER_S:IAER_E) &
            & )

       ! Sedimentation
      !  WRITE(KLUOUT,*) "ACLIB : CALL DUST SEDIM (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( &
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "SEDIM", YPAER_ACLIB%CDUST_SEDIM_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, & 
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PAP=PAP, PAPH=PAPH, PDELP=PDELP, PDZZ=ZDUMMY1, &
            & PTP=PTP, PVERVEL=PVERVEL, PRH=PRH, PRHODREF=PRHODREFM &
            & )

   ENDIF
   DO JLON=KIDIA,KFDIA
    DO JAER=1, YPAER_ACLIB%NAERO
  
        val = PSFAER(JLON,JAER)
  
        ! Check NaN or Inf
        IF (.NOT. IEEE_IS_FINITE(val)) THEN
          PRINT *, "AerMain1 NaN/(AFTER DUST) Inf detected in PSFAER at:", JLON, JAER, val
        END IF
  
        ! Check absurd values (adjust threshold!)
        IF (ABS(val) > 1.E10) THEN
          PRINT *, "(AerMain1 AFTER DUST) Huge value in PSFAER at:", JLON, JAER, val
        END IF
  
      END DO
    END DO
!!!!!!!!!!!!!!!!!!!!!!
!! Sea-salt family  !!
!!!!!!!!!!!!!!!!!!!!!!

   ! Sea-salt family
   IF (YPAER_ACLIB%LSALTON) THEN
       IAER_S = YPAER_ACLIB%NSALT_F0
       IAER_E = YPAER_ACLIB%NSALT_F1
       INAER = IAER_E - IAER_S + 1

       ! Emissions
      !  WRITE(KLUOUT,*) "ACLIB : CALL SALT EMIS (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( & 
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "EMIS", YPAER_ACLIB%CSALT_EMIS_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, &
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PRH=PRH, PRHODREF=PRHODREFM, &
            & PWIND=PWIND, PSST=PSST, PWS1=PWS1, PCI=PCI, &
            & PBARESOIL=PBARESOIL, PLSM=PLSM, PCLK=PCLK, PSNS=PSNS, &
            & PAERDEP=PAERDEP, PAERSCC=PAERSCC, PZ0M=PZ0M, &
            & PSFAER=PSFAER(:,IAER_S:IAER_E) &
            & )

       ! Dry deposition
      !  WRITE(KLUOUT,*) "ACLIB : CALL SALT DRYDEP (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( & 
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "DRYDEP", YPAER_ACLIB%CSALT_DRYDEP_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, &
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PAP=PAP, PAPH=PAPH, PAPHI=PAPHI, PDELP=PDELP, &
            & PRHODREF=PRHODREFM, &
            & PCI=PCI, PLSM=PLSM, PSNS=PSNS, &
            & PSFAER=PSFAER(:,IAER_S:IAER_E) &
            & )

       ! Sedimentation
      !  WRITE(KLUOUT,*) "ACLIB : CALL SALT SEDIM (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( &
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "SEDIM", YPAER_ACLIB%CSALT_SEDIM_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, & 
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PAP=PAP, PAPH=PAPH, PDELP=PDELP, PDZZ=ZDUMMY1, &
            & PTP=PTP, PVERVEL=PVERVEL, PRH=PRH, PRHODREF=PRHODREFM &
            & )

   ENDIF

   DO JLON=KIDIA,KFDIA
    DO JAER=1, YPAER_ACLIB%NAERO
  
        val = PSFAER(JLON,JAER)
  
        ! Check NaN or Inf
        IF (.NOT. IEEE_IS_FINITE(val)) THEN
          PRINT *, "AerMain1 NaN/(AFTER SALT) Inf detected in PSFAER at:", JLON, JAER, val
        END IF
  
        ! Check absurd values (adjust threshold!)
        IF (ABS(val) > 1.E10) THEN
          PRINT *, "(AerMain1 AFTER SALT) Huge value in PSFAER at:", JLON, JAER, val
        END IF
  
      END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Anthropogenic family  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF (YPAER_ACLIB%LANTHON) THEN
       IAER_S = YPAER_ACLIB%NANTH_F0
       IAER_E = YPAER_ACLIB%NANTH_F1
       INAER = IAER_E - IAER_S + 1

       ! Emissions are done in SURFEX
       ! We save here the emission in YAEROUT_ACLIB (diag)

       DO JAER=IAER_S,IAER_E
        !  WRITE(KLUOUT,*) "AER_MAIN1 : EMISSIONS SAVED - ",JAER," - ",YAEROSOLS(JAER)%CAER_NAME
         DO JLON=KIDIA,KFDIA
           YAEROUT_ACLIB%PEMIS(JLON,JAER) = -PSFAER(JLON,JAER)
         ENDDO
       ENDDO

       ! Dry deposition
      !  WRITE(KLUOUT,*) "ACLIB : CALL ANTH DRYDEP (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( & 
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "DRYDEP", YPAER_ACLIB%CANTH_DRYDEP_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, &
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PAP=PAP, PAPH=PAPH, PAPHI=PAPHI, PDELP=PDELP, &
            & PRHODREF=PRHODREFM, &
            & PCI=PCI, PLSM=PLSM, PSNS=PSNS, &
            & PSFAER=PSFAER(:,IAER_S:IAER_E) &
            & )

       ! Sedimentation
      !  WRITE(KLUOUT,*) "ACLIB : CALL ANTH SEDIM (",IAER_S," - ",IAER_E
       CALL AER_ACLIB( &
            & YDCST_ACLIB, YPAER_ACLIB, YAEROSOLS(IAER_S:IAER_E), &
            & "SEDIM", YPAER_ACLIB%CANTH_SEDIM_CODE, &
            & KAERO=INAER, KMODE=0, KCHEM=0, &
            & KAER_START=IAER_S, KIDIA=KIDIA, KFDIA=KFDIA, &
            & KLUOUT=KLUOUT, KVERB=0, &
            & KKA=IDUMMY, KKU=IDUMMY, KKL=-1, & 
            & KLON=KFDIA, KLAT=1, KLEV=KLEV, PTSTEP=PTSTEP, & 
            & PAER=PAER(:,:,IAER_S:IAER_E), &
            & PTENAER=PTENAER(:,:,IAER_S:IAER_E), &
            & YAEROUT_ACLIB=YAEROUT_ACLIB, &
            & PAP=PAP, PAPH=PAPH, PDELP=PDELP, PDZZ=ZDUMMY1, &
            & PTP=PTP, PVERVEL=PVERVEL, PRH=PRH, PRHODREF=PRHODREFM &
            & )

   ENDIF




! Check nan 


DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    DO JAER=1, YPAER_ACLIB%NAERO

       val = PAER(JLON,JLEV,JAER)

       ! Check NaN or Inf
       IF (.NOT. IEEE_IS_FINITE(val)) THEN
         PRINT *, "AerMain1 NaN/Inf detected in PAER at:", JLON,JLEV,JAER, val
       END IF

       ! Check absurd values (adjust threshold!)
       IF (ABS(val) > 1.E10) THEN
         PRINT *, "AerMain1 Huge value in PAER at:", JLON,JLEV,JAER, val
       END IF
       IF (val < 0.) THEN
           PRINT *, "AerMain1 Negative aerosol in PAER at:", JLON,JLEV,JAER,val
         END IF
     END DO
   END DO
 END DO

 DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    DO JAER=1, YPAER_ACLIB%NAERO

       val = PTENAER(JLON,JLEV,JAER)

       ! Check NaN or Inf
       IF (.NOT. IEEE_IS_FINITE(val)) THEN
         PRINT *, "AerMain1 NaN/Inf detected in PTENAER at:", JLON,JLEV,JAER, val
       END IF

       ! Check absurd values (adjust threshold!)
       IF (ABS(val) > 1.E10) THEN
         PRINT *, "AerMain1 Huge value in PTENAER at:", JLON,JLEV,JAER, val
       END IF
     END DO
   END DO
 END DO


DO JLON=KIDIA,KFDIA
  DO JAER=1, YPAER_ACLIB%NAERO

      val = PSFAER(JLON,JAER)

      ! Check NaN or Inf
      IF (.NOT. IEEE_IS_FINITE(val)) THEN
        PRINT *, "AerMain1 NaN/Inf detected in PSFAER at:", JLON, JAER, val
      END IF

      ! Check absurd values (adjust threshold!)
      IF (ABS(val) > 1.E10) THEN
        PRINT *, "AerMain1 Huge value in PSFAER at:", JLON, JAER, val
      END IF

    END DO
  END DO

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_MAIN1',1,ZHOOK_HANDLE)    
!------------------------------------------------------------------------------
END SUBROUTINE AER_MAIN1



! SUBROUTINE check_array_2d(name, arr, dim1, dim2)
!      USE, INTRINSIC :: IEEE_ARITHMETIC
!      USE PARKIND1,             ONLY : JPIM     ,JPRB

!      IMPLICIT NONE
   
!      CHARACTER(LEN=*), INTENT(IN) :: name
!      REAL(KIND=JPRB), INTENT(IN) :: arr(:,:)
!      INTEGER(KIND=JPIM), INTENT(IN) :: dim1, dim2
   
!      INTEGER(KIND=JPIM) :: i, j
!      REAL(KIND=JPRB) :: val
   
!      DO j = 1, dim2
!        DO i = 1, dim1

!          val = arr(i,j)
   
!          IF (.NOT. IEEE_IS_FINITE(val)) THEN
!            PRINT *, "AerMain1 NaN/Inf in ", TRIM(name), " at:", i, j, val
!          END IF
   
!          IF (ABS(val) > 1.E10) THEN
!            PRINT *, "AerMain1 Huge value in ", TRIM(name), " at:", i, j, val
!          END IF

!          IF (val < 0.) THEN
!           PRINT *, "AerMain1 Negative aerosol in ", TRIM(name), " at:", i,j,val
!         END IF
        
!        END DO
!      END DO
   
!    END SUBROUTINE

!    SUBROUTINE check_array_3d(name, arr, dim1, dim2, dim3)
!      USE, INTRINSIC :: IEEE_ARITHMETIC
!      USE PARKIND1,             ONLY : JPIM     ,JPRB

!      IMPLICIT NONE
   
!      CHARACTER(LEN=*), INTENT(IN) :: name
!      INTEGER(KIND=JPIM), INTENT(IN) :: dim1, dim2, dim3
!      REAL(KIND=JPRB), INTENT(IN) :: arr(dim1, dim2, dim3)

!      INTEGER(KIND=JPIM) :: i, j, k
!      REAL(KIND=JPRB) :: val
!      INTEGER(KIND=JPIM) :: nbad
   
!      nbad = 0
   
!      DO JLEV=1,KLEV
!       DO JLON=KIDIA,KFDIA
!         DO k = 1, dim3
!           print*, dim1, dim2, dim3, i, j, k
   
!            val = arr(i,j,k)
   
!            ! Check NaN or Inf
!            IF (.NOT. IEEE_IS_FINITE(val)) THEN
!              PRINT *, "NaN/Inf detected in ", TRIM(name), " at:", i, j, k, val
!              nbad = nbad + 1
!            END IF
   
!            ! Check absurd values (adjust threshold!)
!            IF (ABS(val) > 1.E10) THEN
!              PRINT *, "Huge value in ", TRIM(name), " at:", i, j, k, val
!              nbad = nbad + 1
!            END IF
!            IF (val < 0.) THEN
!                PRINT *, "Negative aerosol in ", TRIM(name), " at:", i,j,k,val
!              END IF
!          END DO
!        END DO
!      END DO
   
!      IF (nbad > 0) THEN
!        PRINT *, "WARNING:", nbad, " problematic values in ", TRIM(name)
!      END IF
   
!    END SUBROUTINE