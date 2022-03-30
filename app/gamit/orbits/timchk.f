      subroutine timchk( jdss,tss,delt,nepoch,jds,ts,jdf,tf,jde,te )

c     Adjust the requested epochs for a T-file to match those available
c     on the NGS SP#1 file

c     R.W. King March 1988
c
      implicit real*8 (a-h,o-z)
c
      integer*4 jds,jdf,jde,jdss, nepoch

c Adjust the start time

      xjdss= dble(jdss) + tss/86400.d0
      xjdst= dble(jds) + ts/86400.d0
      if( xjdst.lt.xjdss ) then
         write(message,'(3a,f12.4,a,f12.4)
     .        'After allowing for interpolation,'
     .       ,' SP3 ephemeris doesn''t start early enough'
     .       ,'--Requested JD=',jdss,'  Modified SP3 JD=',xjdst
         call report_stat('FATAL','ORBITS','timchk',' ',message,0)
      endif
      dtab= dint( (xjdst-xjdss)*86400.d0/delt )
      dt  = dtab*delt
      jds = jdss
      ts= tss
      call timinc( jds,ts,dt )
CD     write(6,11) jds,ts,jdss,tss,dtab,dt
   11 format(1x,'In TIMCHK, JDS,TS,JDSS,TSS,DTAB,DT=',/
     1      ,1x,i7,f13.6,1x,i7,f13.6,f7.2,f15.6)

c  Adjust the stop time
      xjdfs= xjdss + delt*dble(nepoch-1)/86400.d0
      xjdft= dble(jdf) + tf/86400.d0
      if( xjdft.gt.xjdfs )then
        write(message,'(3a,f12.4,a,f12.4)
     .        'After allowing for interpolation,'
     .       ,' SP3 ephemeris doesn''t go late enough'
     .       ,'--Requested JD=',jdss,'  Modified SP3 JD=',xjdst
         call report_stat('FATAL','ORBITS','timchk',' ',message,0)


        write(6,20) xjdft,xjdfs
   20   format(1x,'NGS ephemeris doesn''t go late enough,',/
     1        ,4x,'Requested JD=',f12.4,'  NGS JD=',f12.4)
        stop
      endif
      dtab= dint( (xjdft-xjdss)*86400.d0/delt )
      dt= dtab*delt
      jdf= jdss
      tf = tss
      call timinc( jdf,tf,dt )
CD     write(6,21) jdf,tf,dt
   21 format(1x,'In TIMCHK, JDF,TF,DT=',i7,f13.6,f15.6)


c Adjust the ICs epoch

      xjde= dble(jde) + te/86400.d0
      dtab= dint( (xjde-xjdss)*86400.d0/delt )
      dt  = dtab*delt
      jde = jdss
      te  = tss
      call timinc( jde,te,dt )
CD     write(6,31) jde,te,jdss,tss,dtab,dt
   31 format(1x,'In TIMCHK, JDE,TE,JDSS,TSS,DTAB,DT=',/
     1      ,1x,i7,f13.6,1x,i7,f13.6,f7.2,f15.6)


      return
      end
