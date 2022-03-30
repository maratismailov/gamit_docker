      real*8 function clkdif
     .       (jdoy1,jhour1,jmin1,sec1,
     .        jdoy0,jhour0,jmin0,sec0)
      real*8 sec0,sec1
      integer jdoy1,jdoy0
      integer jhour1,jhour0
      integer jmin0,jmin1

c     given two times t1 and t0 in jdoy,jhour,jmin,sec
c     return t1 - t0 in seconds

      clkdif = dble(jdoy1  - jdoy0)   * 86400.0d0
     .       + dble(jhour1 - jhour0)  *  3600.0d0
     .       + dble(jmin1  - jmin0)   *    60.0d0
     .       + sec1 - sec0

      return

      end
