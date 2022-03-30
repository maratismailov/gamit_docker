      Integer*8 Function itimdif ( itime1, itime2 )
c
c     Return difference in seconds between epoch 1 minus epoch 2
c     handles going over days and positive and negative diffs.  If the
c     difference exceeds integer*4 capacity (2**30 sec = 12000 days = 34 years),
c     set the difference at 30 years (=11000 days = ).  Also trap open-ended 
c     dates (from station.info) designated by 9999 999

c     R. King  7 January 2003  

c     Input times are integer year, day-of-year, hour, minutes, seconds
             
      implicit none

      integer*4 itime1(5),itime2(5),jd1,jd2,sod1,sod2,id,im,spd,julday
     .        , iyr,idoy

C     seconds per day
      data        spd/86400/
         
      iyr = itime1(1)
      idoy = itime1(2)
      if( iyr.eq.9999 ) then   
        iyr = 2100
        idoy = 0    
      endif
      call monday(idoy,im,id,iyr) 
cd  cc  print *,'ITIMDIF 1 iyr idoy im id ',iyr,idoy,im,id 
      jd1 = julday(im,id,iyr) 
      sod1 = 3600*itime1(3) + 60*itime1(4) + itime1(5)   
      iyr = itime2(1)
      idoy = itime2(2)
      if( iyr.eq.9999 ) then   
        iyr = 2100
        idoy = 0 
      endif    
      call monday(idoy,im,id,iyr)    
cd      print *,'ITIMDIF 2 iyr idoy im id ',iyr,idoy,im,id 
      jd2 = julday(im,id,iyr) 
      sod2 = 3600*itime2(3) + 60*itime2(4) + itime2(5) 
      if( abs(jd1-jd2).gt.11000 ) then
        itimdif = sign(11000,(jd1-jd2))*spd
      else
        itimdif = spd*(jd1-jd2) + sod1-sod2 
      endif
cd      print *,'ITIMDIF jd1 jd2 sod1 sod2 '
cd     .     , itimdif,jd1,jd2,sod1,sod2 
      return
      end
