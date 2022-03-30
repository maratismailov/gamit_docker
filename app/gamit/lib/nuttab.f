      SUBROUTINE NUTTAB( INUT,PJD,OBLQ,DPSI,DEPS,RNUT )
C
C     NUTTAB takes JD as an argument and return the nutations
C     in longitude and obliquity (in radians) and the nutation matrix.
C     Values of DEPS AND DPSI are interpolated from tables.
C     R.I. Abbot - June 1984
C     Merged orbits/model version moved to libary by R. King - January 1992
c     Sense of matrix (from /orbits) opposite from original /model version
c     Changed to call MHB_2000 for nutations if the nbody file  exists
c     (for the lunar and solar ephemerides); otherwise use the old nutabl.
c     Note that MBH_2000 requires true JD not PEP JD and returns milliarcseconds
c     not arcsecons.  R. King 180328
C
C     Input
C
C        INUT   Unit number of nutation table
C        PJD    Julian date in TDT 
C        OBLQ   Mean obliquity of date
C
C      Output
C
C        DPSI   Delta psi
C        DEPS   Delta epsilon
C        RNUT   Nutation matrix

      implicit none

      integer*4 inut                                          

      logical fcheck

      REAL*8 RNUT(3,3),RNA(3,3),RNB(3,3),RNC(3,3),WORK(3,3)
     .     , pc,ps,pi,casr,deps,dpsi,oblq,pjd,tjd,cobliq,sobliq

c     Arguments for MHB_2000 not used
      real*8 dpsi_ls,deps_ls,dpsi_plan,deps_plan,dpsi_fcn,deps_fcn
     .     ,dpsi_prec,deps_prec

      PI= 4.D0*DATAN(1.D0)
      CASR= PI/180.D0/3600.D0

c        Obtain nutation parameters either from the nuttabl. or by
c        evalutating the formulas directly. 
      
      if( fcheck('nbody') ) then  
        tjd = pjd - 0.5d0 
        call MHB_2000(tjd,dpsi_ls,deps_ls
     .                   ,dpsi_plan, deps_plan
     .                   ,dpsi_fcn,  deps_fcn
     .                   ,dpsi_prec, deps_prec
     .                   ,dpsi      ,deps )
        dpsi = dpsi/1000.d0
        deps = deps/1000.d0 
cd        print *,'from MHB dpsi deps ',dpsi,deps
      else  
        CALL NUTRED( INUT,pjd,DPSI,DEPS )      
cd        print *,'from NUTRED dpsi deps ',dpsi,deps 
      endif

      DPSI=DPSI*CASR
      DEPS=DEPS*CASR
C
C
C        Auxiliary quantities for the computation of the nutation matrix

      COBLIQ=DCOS(OBLQ)
      SOBLIQ=DSIN(OBLQ)
      PC=COBLIQ*DPSI
      PS=SOBLIQ*DPSI


c      Calculation of the nutation matrix
C
c     RNUT(1,1) = 1.D0 - 0.5D0*DPSI**2
c     RNUT(2,1) = -PC
c     RNUT(3,1) = -PS
c     RNUT(1,2) = PC - DEPS*PS
c     RNUT(2,2) = 1.D0 - 0.5D0*(DEPS**2 + PC**2)
c     RNUT(3,2) = -DEPS - 0.5D0*PC*PS
c     RNUT(1,3) = PS + DEPS*PC
c     RNUT(2,3) = DEPS - 0.5D0*PC*PS
c     RNUT(3,3) = 1.D0 - 0.5D0*(DEPS**2 + PS**2)


c       Calculate the nutation matrix (Reference:  Mueller, p. 75)

      Call ROTMAT(OBLQ,1,RNA)
      Call ROTMAT(-DPSI,3,RNB)
      Call ROTMAT(-OBLQ-DEPS,1,RNC)
      Call MATMPY(RNB,RNA,WORK,3,3,3)
      Call MATMPY(RNC,WORK,RNUT,3,3,3)

      RETURN
      END
