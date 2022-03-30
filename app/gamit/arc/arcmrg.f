Copyright (c) Massachusetts Institute of Technology and The University of
California, 1994.  All rights reserved.
      subroutine arcmrg

c     merge individual satellite ephemerides into one
c     Rick Abbot - November 1984

      implicit none

      include '../includes/dimpar.h'
      include '../includes/units.h'
      include '../includes/arc.h'

      character*80 head,comt(3)
      character*2  buf2
      character*1  upperc

      real*8 tee,t00,tff,satics2,tabs2,dpvdm2

      integer*4 ite,it0,itf,nepoch,nintrp,nsat1
     .        , iunit,nprtl,ncrd,ioerr,isat,iepoch,i,j

      dimension tee(3),t00(3),tff(3)
      dimension ite(3),it0(3),itf(3)
      character*4 sname4(4,maxsat)
      dimension satics2(maxorb,maxsat),tabs2(6,maxsat)
     1        , dpvdm2(maxyt2,maxsat)
                  

c     Read individual integrations from scratch units 40-(40+nsats)
c     Write the T-file onto unit 18
                  
c------ Read the headers to get common plus each satellite name    

      nprtl = 0  
      do isat=1,nsats
         iunit=isat+39
         read(iunit) head,ite,tee,it0,t00,itf,tff,delt,nepoch
     .             , ut1utc,xpole,ypole     
         write  (buf2,'(a2)') head(79:80)
         read (buf2,'(i2)') nics
         read (iunit) comt,nsat1,nintrp
     .           , (sname4(i,isat),i=1,4),(satics2(i,isat),i=1,nics)
         if( upperc(apar).eq."V" ) then
           ncrd = 6
           nprtl = nintrp - 6
         else      
           ncrd = 3
           nprtl = nintrp - 3
         endif
       enddo
          

c------- Write the T-file header
                               
cd      print *,'ARCMRG iut nsat1 nsats ',iut,nsat1,nsats
      write (iut,err=993,iostat=ioerr)
     .   head,ite,tee,it0,t00,itf,tff,delt,nepoch,ut1utc,xpole,ypole
      write (iut,err=993,iostat=ioerr)
     .   comt,nsats,nintrp,
     . ((sname4(i,isat),i=1,4),(satics2(i,isat),i=1,nics),isat=1,nsats)
            

c------- Read and write the data record for each epoch
             
cd      print *,'ARCMRG nepoch ',nepoch           
      do iepoch=1,nepoch
         do isat=1,nsats
            iunit=isat+39                  
cd            print *,'iepoch isat iunit ',iepoch,isat,iunit
            if( nprtl.eq.0 ) then 
               read (iunit,end=997,err=995,iostat=ioerr)
     .              (tabs2(i,isat),i=1,ncrd)     
            else
               read (iunit,end=997,err=995,iostat=ioerr)
     .               (tabs2(i,isat),i=1,ncrd)
     .             , (dpvdm2(i,isat),i=1,nprtl)
            endif 
         enddo
         if( nprtl.eq.0 ) then
            write(iut,err=993,iostat=ioerr)
     .           ((tabs2(i,j),i=1,ncrd),j=1,nsats)
         else
            write (iut,err=993,iostat=ioerr)
     .            ( (tabs2(i,j),i=1,ncrd),(dpvdm2(i,j),i=1,nprtl)
     .             , j=1,nsats )
         endif  
      enddo

c------- Close the scratch units

      do isat=1,nsats
         iunit=isat+39
         close (unit=iunit,status='delete')
      enddo

c------- Exit
      goto 999   
                            

c------- Errors encountered 

  993 continue
      call report_stat('FATAL','ARC','arcmrg',' ',
     .'Error writing T-file',ioerr)
  995 continue
      call report_stat('FATAL','ARC','arcmrg',' ',
     .'Error reading scratch file',ioerr)
  997 continue
      call report_stat('FATAL','ARC','arcmrg',' ',
     .'Error premature end of scratch file',ioerr)

  999 return

      end
