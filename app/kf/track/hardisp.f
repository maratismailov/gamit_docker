      SUBROUTINE HARDISP(yr_o,mon_o,dy_o,h_o,m_o,s_o,
     .                   in_unit,dNEU_otl)
* MOD AZ 190305: subroutine called in theory.f
* Input: date_otl(1) date_otl(2) date_otl(3) date_otl(4)
*        date_otl(5) seconds_otl_int blq_unit
* Output: dNEU_otl
* Original comments of IERS routines kept for further elaboration
*
!*      PROGRAM HARDISP
!+
!
!  - - - - - - - - - - -
!   H A R D I S P
!  - - - - - - - - - - -
!
!  This program is part of the International Earth Rotation and
!  Reference Systems Service (IERS) Conventions software collection.
!
!  This program reads in a file of station displacements in the BLQ
!  format used by Scherneck and Bos for ocean loading, and outputs a
!  time series of computed tidal displacements, using an expanded set
!  of tidal constituents, whose amplitudes and phases are found by
!  spline interpolation of the tidal admittance.  A total of 342
!  constituent tides are included, which gives a precision of about
!  0.1%.
!
!  In general, Class 1, 2, and 3 models represent physical effects that
!  act on geodetic parameters while canonical models provide lower-level
!  representations or basic computations that are used by Class 1, 2, or
!  3 models.
! 
!  Status:  Class 1 model
!
!     Class 1 models are those recommended to be used a priori in the
!     reduction of raw space geodetic data in order to determine
!     geodetic parameter estimates.
!     Class 2 models are those that eliminate an observational
!     singularity and are purely conventional in nature.
!     Class 3 models are those that are not required as either Class
!     1 or 2.
!     Canonical models are accepted as is and cannot be classified as
!     a Class 1, 2, or 3 model.
!
!  Given:
!     User provided input ocean loading coefficients (Note 1)
!
!  Returned:
!     DU         d      Radial tidal ocean loading displacement (Note 2)
!     DW         d      West tidal ocean loading displacement (Note 2)
!     DS         d      South tidal ocean loading displacement (Note 2)
!
!     :------------------------------------------:

      IMPLICIT NONE
* MOD AZ 190305: include Header for shared variable "usr_interval"
      include 'track_com.h'

      INTEGER yr_o, mon_o, dy_o, h_o, m_o, s_o
      INTEGER num_o,in_unit
      real*8 dNEU_otl(3)
      INTEGER*4 I,IDAY,IDT_otl,IMONTH,IRNT,IRHI,IRLI,
     .        IT,NB,NL,NP,NT,NTIN,
     .        KK,NTOUT,MDAY
      PARAMETER (NL=600)
      PARAMETER (NT=342)
      PARAMETER (NTIN=11)

      CHARACTER*40 DUMM
      REAL AMP,AS,AW,AZ,DS,DW,DZ,HCS,HCW,HCZ,PHASE,TAMP,TPH,WF,SAMP
      DOUBLE PRECISION F,PZ,PS,PW,SCR
      DOUBLE PRECISION DR,PI

      DIMENSION TAMP(3,NTIN),TPH(3,NTIN)
      DIMENSION IDT_otl(6,NTIN),AMP(NTIN),PHASE(NTIN)
      DIMENSION AZ(NT),PZ(NT),HCZ(2*NT)
      DIMENSION AS(NT),PS(NT),HCS(2*NT)
      DIMENSION AW(NT),PW(NT),HCW(2*NT)
      DIMENSION DZ(NL),DS(NL),DW(NL)
      DIMENSION F(NT),SCR(3*NT),WF(NT)
      COMMON/DATE_OTL/IT(5)
      DATA DR/0.01745329252D0/,IRLI/1/
      PARAMETER ( PI = 3.1415926535897932384626433D0 ) 


      DATA IDT_otl/
     .  2, 0, 0, 0, 0, 0,   2, 2,-2, 0, 0, 0,   2,-1, 0, 1, 0, 0,
     .  2, 2, 0, 0, 0, 0,   1, 1, 0, 0, 0, 0,   1,-1, 0, 0, 0, 0,
     .  1, 1,-2, 0, 0, 0,   1,-2, 0, 1, 0, 0,   0, 2, 0, 0, 0, 0,
     .  0, 1, 0,-1, 0, 0,   0, 0, 2, 0, 0, 0/

      IT(1)=yr_o
      IT(2) = dy_o + MDAY(IT(1),mon_o)
      IT(3) = h_o
      IT(4) = m_o
      IT(5) = s_o
* MOD AZ 190305: IRNT defined as 1 since epoch-wise OTL needed
      IRNT = 1
      SAMP = usr_interval
*      write(*,*) IT(1),IT(2),IT(3),IT(4),IT(5),SAMP
*      write(*,*) 'Okay before BLQ read-in' 
      DO I=1,3
        READ(in_unit,108) (TAMP(I,KK),KK=1,NTIN)
 108    FORMAT(1X,11F7.5)
      ENDDO
      DO I=1,3
        READ(in_unit,110) (TPH(I,KK),KK=1,NTIN)
 110    FORMAT(1X,11F7.1)

        DO KK=1,NTIN
          TPH(I,KK)=-TPH(I,KK)
        ENDDO
      ENDDO
      
*      write(*,*) 'Okay till after BLQ read-in'
*      write(*,*) TAMP(3,11)
*      write(*,*) TPH(3,11)
     
      DO I=1,NTIN
        AMP(I)=TAMP(1,I)
        PHASE(I)=TPH(1,I)
      ENDDO
      CALL ADMINT(AMP,IDT_otl,PHASE,AZ,F,PZ,NTIN,NTOUT)
      DO I=1,NTIN
        AMP(I)=TAMP(2,I)
        PHASE(I)=TPH(2,I)
      ENDDO
      CALL ADMINT(AMP,IDT_otl,PHASE,AW,F,PW,NTIN,NTOUT)
      DO I=1,NTIN
        AMP(I)=TAMP(3,I)
        PHASE(I)=TPH(3,I)
      ENDDO
      CALL ADMINT(AMP,IDT_otl,PHASE,AS,F,PS,NTIN,NTOUT)

*      write(*,*) 'Okay till after call ADMINT'
      
      DO I=1,NTOUT
        PZ(I) = DR*PZ(I)
        PS(I) = DR*PS(I)
        PW(I) = DR*PW(I)
        F(I) = SAMP*PI*F(I)/43200.D0
        WF(I) = F(I)
      ENDDO


 11   IRHI = MIN(IRLI+NL-1,IRNT)
      NP = IRHI - IRLI + 1

      DO I=1,NT
        HCZ(2*I-1) = AZ(I)*DCOS(PZ(I))
        HCZ(2*I)  = -AZ(I)*DSIN(PZ(I))
        HCS(2*I-1) = AS(I)*DCOS(PS(I))
        HCS(2*I)  = -AS(I)*DSIN(PS(I))
        HCW(2*I-1) = AW(I)*DCOS(PW(I))
        HCW(2*I)  = -AW(I)*DSIN(PW(I))
      ENDDO
      
*      write(*,*) 'Okay till before call RECURS'
      
      CALL RECURS(DZ,NP,HCZ,NTOUT,WF,SCR)
      CALL RECURS(DS,NP,HCS,NTOUT,WF,SCR)
      CALL RECURS(DW,NP,HCW,NTOUT,WF,SCR)
* MOD AZ 190305: Output returned in NEU
      dNEU_otl(1) = -DS(1)
      dNEU_otl(2) = -DW(1)
      dNEU_otl(3) = DZ(1)

*      write(*,*) 'Okay till after call RECURS'
*      write(*,*) dNEU_otl(1),' ',dNEU_otl(2),' ',dNEU_otl(3) 
      IF(IRHI.EQ.IRNT) RETURN
      IRLI = IRHI + 1

      DO I=1,NT
        PZ(I) = DMOD(PZ(I) + NP*F(I),2.D0*PI)
        PS(I) = DMOD(PS(I) + NP*F(I),2.D0*PI)
        PW(I) = DMOD(PW(I) + NP*F(I),2.D0*PI)
      ENDDO
      GO TO 11


      END
