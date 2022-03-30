      Subroutine even_minute(yr,doy,sod,nyr,ndoy,nsod )

c     Get the correct year and day for calls to session.info and station.info 
c     by rounding the (start) time to the even minute.  This mainly serves to 
c     get the day right in the case of X- or C-files in which the sampling is 
c     several seconds prior to the even minute (e.g., TI firmware sampling at 
c     59.00 or 59.08 GPST and MiniMac and Rogue firmware sampling on even UTC 
c     epochs.  

c     R. King 6 May 1997

c       Input
c             yr      I*4   actual year of epoch
c             doy     I*4   actual day of epoch
c             sod     R*8   seconds-of-day of epoch

c       Output is input rounded to the nearest minute
c             nyr     I*4   year 
c             ndoy    I*4   day of year
c             nsod    I*4   seconds-of-day 

      implicit none

      integer*4 yr,doy,nyr,ndoy,nsod,month,day,jd,julday,min

      real*8 sod
           
      call monday(doy,month,day,yr)
      jd = julday(month,day,yr)  
      min = nint(sod/60.) 
      if( min.ge.1440 ) then
         min = 0
         jd = jd + 1
      endif
      nsod = min*60
      call dayjul( jd,nyr,ndoy )
        
      return
      end        
      



