MODULE M_TACTIC_CGROWTH

    CONTAINS

    SUBROUTINE TACTIC_CGROWTH &
        !---input
        & ( KIDIA, KFDIA, KLON, KLEV, &
        &   PXRATE, &
        &   PAERPHO , PITAERPHO, PTSPHY  ,&
        !---output
        &   PTAERPHI, PTAERPHO )
    
        !*** * TACTIC_CGROWTH* - GROWTH OF O.M. AND B.C. AEROSOLS
        
        !**   INTERFACE.
        !     ----------
        !   *TACTIC_CGROWTH* IS CALLED FROM *AER_PHY3*.
        
        !     AUTHOR.
        !     -------
        !        JEAN-JACQUES MORCRETTE  *ECMWF*
        !        FROM O.BOUCHER's cgrowth 
        
        !     SOURCE.
        !     -------
        
        !     MODIFICATIONS.
        !     --------------
        !        ORIGINAL : 2006-11-21
        !-----------------------------------------------------------------------
        
        USE PARKIND1  ,ONLY : JPIM     ,JPRB
        USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
        
        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        
        !*       0.1   ARGUMENTS
        !              ---------
        
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
        INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
        
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PXRATE
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
        REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERPHO(KLON,KLEV) , PITAERPHO(KLON,KLEV)
        REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERPHI(KLON,KLEV), PTAERPHO(KLON,KLEV)
        
        
        !*       0.5   LOCAL VARIABLES
        !              ---------------
        
        INTEGER(KIND=JPIM) :: JL, JK
        REAL(KIND=JPRB) :: ZAERCONV
        
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_CGROWTH:TACTIC_CGROWTH',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
        
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZAERCONV = (PAERPHO(JL,JK) + PTSPHY * PITAERPHO(JL,JK)) * PXRATE

        ! Limit by available reservoir
            ZAERCONV = MIN(ZAERCONV, PAERPHO(JL,JK))
            
            PTAERPHI(JL,JK) =  ZAERCONV
            PTAERPHO(JL,JK) = -ZAERCONV    
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M_TACTIC_CGROWTH:TACTIC_CGROWTH',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
        END SUBROUTINE TACTIC_CGROWTH

END MODULE M_TACTIC_CGROWTH
