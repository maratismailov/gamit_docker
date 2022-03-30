CTITLE    ................................................................
 
      SUBROUTINE epoc_8(IM,ID,IY,IHR,IMIN,TIM)
C
c mod tah 850617 Changed to use real*8 for use in SOLVK
c
C     **************************************************************
C     *                                                            *
C     *     D. ROBERTSON             AUG. 1976                     *
C     *                                                            *
C     *  EPOC IS A SUBROUTINE WHICH WILL CONVERT A JULIAN DATE     *
C     *  INTO MONTH, DAY, YEAR, HOUR, AND MINUTE.                  *
C     *                                                            *
C     **************************************************************
C
      real*8 TIM,FTIM,FRAC
 
      integer*4 im, id, iy, ihr, imin, ifct, it, itime
 
C
      IFCT = TIM/1.0D4
      FTIM = TIM - IFCT * 1.0D4
C
      IT = FTIM - 0.5D0
      FRAC = FTIM - 0.5D0 - IT
      FTIM = FTIM - FRAC + IFCT * 1.0D4
C
      CALL MDYJL_8(IM,ID,IY,ITIME,FTIM)
C
      IHR = FRAC*24.0D0
      FRAC = FRAC - IHR / 24.0D0
      IMIN = FRAC*1440.0D0
C
      RETURN
      END
 
