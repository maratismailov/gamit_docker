      subroutine writee( lu,icall,nhead,head,trans_sow
     .                 , nprn,iewkn,bephem,bclock,subfr1 )

C     Write out a version 2 RINEX navigation (GAMIT 'e') file on unit lu (assumed open)
c
c     Input 
c       lu      : logical unit number (assumed open)
c       icall   : Called to write header(0) or data (1) records
c       nhead     : number of header records
c       trans_sow : Transmission time seconds of GPS week 
c       nprn     : Satellite ID number (PRN)
c       iewkn    : Week number at reference epoch
c       bephem   : Broadcast ephemeris elements 
c       bclock   : Broadcast clock parameters
c       subfr1   : Subframe 1 parameters

C     EPHEMERIS PARAMETERS bephem(i)
C      1   xetoe
C      2   xem0
C      3   xedn
C      4   xeart     square root of semi-major axis, in m**1/2
C      5   xeecc
C      6   xei0
C      7   xeidt
C      8   xeom0
C      9   xeomd
C     10   xew
C     11   xecuc
C     12   xecus
C     13   xecrc
C     14   xecrs
C     15   xecic
C     16   xecis
c
c     CLOCK PARAMETERS bclock(i)
c      1   sow
c      2   xeaf0
c      3   xeaf1
c      4   xeaf2
c
c     SUBFRAME 1 PARAMETERS subfr(i)
c      1     L2 code flag cflgl2
c      2     SV accuracy
c      3     SV health
c      4     Age of Data Clock
c      5     L2 P data flag
c      6     Group delay differential
c      7     Age of Data (ephemeris)

c     TRANSMISSION TIME trans_sow in RINEX2 line 7 derived from 
c       Z-count in Hand Over Word (HOW)

      implicit none
       
      character*80 prog_name

C     Number of elements in array

      real*8 bephem(16),bclock(6),subfr1(8)

      integer*4 it,lu,iday,jdoy,nprn,imonth,ihr
     .        , min,iyear,iyr2,iewkn,icall,i
     .        ,irnxver,len,rcpar
                   
c     RINEX navigation file variables
      real*8
     .         sow,sec,utcoff
     .       , xeaf0,xeaf1,xeaf2,aode,xecrs,xedn,xem0,xecuc,xeecc
     .       , xecus,xeart,xetoe,xecic,xeom0,xecis,xei0,xecrc,xew,xeomd
     .       , xeidt,cflgl2,weekno,pflgl2,svaccr,svhlth,tgd,aodc
     .       , trans_sow,spare1,spare2,spare3     

c     local variables
      integer maxhed,nhead 
      parameter(maxhed=20)  
      character*80 head(maxhed)
      
c     set parameter values for RINEX 2                 
      parameter(irnxver=2)
      parameter(spare1=0.d0,spare2=0.d0,spare3=0.d0)    

c     get program name calling writee
      len = rcpar(0,prog_name)         
      if( prog_name(1:2).eq.'BC' .or. prog_name(1:2).eq.'bc' )
     .    prog_name = 'ORBITS'  
               

c     write the header lines

      if (icall.eq. 0) then 
                  
         do i = 1,nhead
            write(lu,'(a80)') head(i)
         enddo


c      write the message lines

       else if (icall.eq.1 ) then

c        Line 1 epoch and clock  parameters  
         sow   =  bclock(1)
         it = 4
         call timcon (it,iewkn,sow,iyear,jdoy,ihr,min,sec,utcoff) 
         call monday(jdoy,imonth,iday,iyear)
         sow   =  bclock(1)  
         xeaf0 =  bclock(2)  
         xeaf1 =  bclock(3)  
         xeaf2 =  bclock(4) 
         iyr2 = mod(iyear,100) 
         write(lu,'(i2,5i3,f5.1,3d19.12)')  
     .              nprn,iyr2,imonth,iday,ihr,min,sec
     .            , xeaf0,xeaf1,xeaf2

c        Ephemeris paramters   
c        Line 2
         aode    =   subfr1(7 ) 
         xecrs   =   bephem(14) 
         xedn    =   bephem( 3) 
         xem0    =   bephem( 2) 
         write(lu,20) aode,xecrs,xedn,xem0  
          
   20    format(3x,4d19.12)
               
c        Line 3
         xecuc   =   bephem(11) 
         xeecc   =   bephem( 5) 
         xecus   =   bephem(12) 
         xeart   =   bephem( 4) 
         write(lu,20) xecuc,xeecc,xecus,xeart


c        Line 4
         xetoe   =   bephem( 1) 
         xecic   =   bephem(15)  
         xeom0   =   bephem( 8)  
         xecis   =   bephem(16)  
         write(lu,20)  xetoe,xecic,xeom0,xecis
         
c        Line 5
         xei0    =   bephem( 6) 
         xecrc   =   bephem(13) 
         xew     =   bephem(10) 
         xeomd   =   bephem( 9) 
         write(lu,20)  xei0,xecrc,xew,xeomd

c        Line 6
         xeidt   =   bephem( 7) 
         cflgl2  =   subfr1(1)  
         weekno  =   iewkn      
         pflgl2  =   subfr1(5)   
         write(lu,20)  xeidt,cflgl2, weekno,pflgl2

c        Line 7
         svaccr  =   subfr1(2)  
          svhlth  =   subfr1(3)  
         tgd     =   subfr1(6)  
         aodc    =   subfr1(4)    
         write(lu,20)  svaccr,svhlth,tgd,aodc
          
c        Line 8  Time of transmission and spares 
         write(lu,20) trans_sow,spare1,spare2,spare3

      else
         call report_stat('FATAL',prog_name,'writee',' '
     .            ,'Bad icall value ',0)

      endif

      return
      end
