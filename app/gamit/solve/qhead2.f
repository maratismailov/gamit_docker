c
      Subroutine QHEAD2( phi,phi2,mode )
c
c     mode = 1: output effective phase data on screen
c     mode = 2: print double difference observations by satellite
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'models.h'

      character*80 afmt

      integer mode,istat,jsat,jj

      real*8 phi(maxsit,maxsat),phi2(maxsit,maxsat)
                 

      if (mode.eq.1) then
         if( logprt ) then
           if(l2flag.eq.0.or.l2flag.eq.-1) write(6,20)
           if(l2flag.eq.-2) write(6,40)
           write(afmt,30) nsat
           write(6,afmt)((phi(istat,jsat),jsat=1,nsat),istat=1,nsite)
           if(l2flag.gt.0) then
              write (6,40)
              write(6,afmt)
     .        ((phi2(istat,jsat),jsat=1,nsat),istat=1,nsite)
           endif   
         endif
         goto 100
      endif
 20   format(/,'  l1 phase residuals - station/satellite',/)
 30   format ('(',1x,i2,'f12.3',')')
 40   format(/,'  l2 phase residuals - station/satellite',/)
c
      if (mode.eq.2) then
         write (10,*) ' '
         if( logprt ) write ( 6,*) ' '
         write (10,45)
         if( logprt ) write ( 6,45)
 45      format('           C-file       Elev      Number of double '
     .         ,'differences for each satellite PRN ')
         write (10,50) (isprn(jj),jj=1,nsat)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
         if( logprt ) write ( 6,50) (isprn(jj),jj=1,nsat)
 50      format (23x,'Cutoff',50i5)
         do 60 istat=1,nsite
            write(10,70) istat,cfiln(istat),elevcut(istat)
     .                  ,(iuse(istat,jj),jj=1,nsat)
            if( logprt ) write( 6,70) istat,cfiln(istat),elevcut(istat)
     .                  ,(iuse(istat,jj),jj=1,nsat)
 60      continue
         write (10,*) ' '
         if( logprt ) write ( 6,*) ' '
      endif
 70   format(1x,'OBS',1x,i3,1x,a14,1x,f5.2,50i5)
c
 100  continue
c
      return
      end


