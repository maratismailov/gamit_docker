c     Program to create sig_new entries for GLOBK from RWK notes   
c     This version works for long-term and downweight a range of years

c     The name is short for quick typing    
c     The command is 'lrw site yy1 yy2 hor ver'
c         where yy can be 2 or 4 digits, the dates are start and stop, and the last two
c         entrires are horizontal and vertical noise in mm    
c     Mod rwk 090903: if the yy entry is 'all', then no dates entered and next entry is 'hor'
c     The output file (tempout) is appended   

c     run

      implicit none
                     
      integer iarg,iclarg,yr1,yr2,mo1,mo2,day1,day2
     .      , hr1,hr2,mn1,mn2,ioerr
                     
      character*2 amo,aday
      character*3 adoy
      character*4 ayr
      character*8 site,ahsig,avsig 

      real*4 hsig,vsig                                  
      
      iarg = iclarg(1,site)

      if(iarg.eq.0) then
        write(*,'(a)') 
     .    'Runstring: lrw [site] [yr-start] [yr-end] [Hor] [Vert]' 
        write(*,'(a)') 'Enter sigmas in mm '
      else      
        open(unit=1,file='tempout',status='unknown',access='append'
     .               ,iostat=ioerr)    
        if( ioerr.ne.0 ) then
          write(*,'(a,i5)') 'Error opening output file ',ioerr
           stop
        endif
        iarg = iclarg(2,ayr)   
        read(ayr,*) yr1
        call fix_y2k(yr1)
        iarg = iclarg(3,ayr)    
        read(ayr,*) yr2    
        call fix_y2k(yr2)
        iarg = iclarg(4,ahsig)   
        read(ahsig,*) hsig
        iarg = iclarg(5,avsig)  
        if( iarg.ne.0 ) then
            read(avsig,*) vsig
        else
          vsig = 0. 
        endif
        mo1 = 1
        day1 = 1
        mo2 = 12
        day2 = 31
        hr1 = 0
        mn1 = 0
        hr2 = 24
        mn2 = 0  
        write(1,'(1x,a,1x,a8,3f7.3,2(1x,i5,4i3))') 
     .        'sig_neu'
     .      , site,hsig/1000.,hsig/1000.,vsig/1000.
     .      , yr1,mo1,day1,hr1,mn1,yr2,mo2,day2,hr2,mn2 
      endif
      stop
      end

c************************************************************************8

      integer*4 function iclarg ( iel, arg )

c     return iel'th argument from the command line in arg
c     length of the string is returned in iclarg

c     This is the Sun version

c     equivalent to rcpar of HP1000 library

      integer*4 iel, rcpar
      character*(*) arg

c     call getarg( iel, arg )
c     iclarg = nblen( arg )
* MOD TAH 010610: Use rcpar which has been modfied to
*     account for different starting argument numbers
      iclarg = rcpar( iel, arg)
      
      return
      end   

c******************************************************************

CTITLE RCPAR
 
      integer*4 function rcpar( iel, arg )

 
*     Routine to emulate RCPAR using the igetarg UNIX subroutine
*     Modified to use getarg
 
*         iel       - Element of runstring to get
*       igetarg     - UNIX routine to read runstring
*       len_arg     - Length of arg.
*       trimlen     - Get length of string
*       offset      - Offset to be applied to the passsed element
*                    (0 is assumed to program name, 1 first argument)
 
      integer*4 iel, len_arg, trimlen, offset
 
*             arg   - Arg of runstring
 
      character*(*) arg
      character*4 test_arg
      
      data offset / -1 /
 
****  Get length of argument and runstring
* MOD TAH 010610: To see where the count starts for getarg
      if( offset.lt.0 ) then
          call getarg(0, test_arg)
	  len_arg = trimlen(test_arg)
	  if( len_arg.eq.0 ) then
	      offset = 1
	  else
	      offset = 0
	  end if
      end if
      
      len_arg = LEN(arg)
      call getarg( iel+offset, arg )
      rcpar = trimlen( arg )
 
***** Thats all
      return
      end
 
c*******************************************************************

CTITLE 'TRIMLEN'
 
 
      integer*4 function trimlen(string)
 

* MOD TAH 050129: Adding counting null's as well as blanks at the 
*     end of strings.  Useful for strings returned from C that have
*     nulls at the end.
 
*     Routine to return the length of the used portion of string.
*     Length is with trailing blanks removed.
 
*         len_string    - declared length of string
*         i             - Loop counter
 
      integer*4 len_string, i
 
*       blanks          - Indicates blanks are still being
*                       - found
 
      logical blanks
 
*             string    - the string whose length we want
 
      character*(*) string
 
***** Get full length of string, and work backwards
 
      len_string = LEN(string)
      i = len_string
 
      blanks = .true.
 
*     Scan string from end to front
      do while ( blanks )
          if( i.gt.0 ) then
* MOD TAH 050129: Add or check on null
              if( string(i:i).eq.' ' .or.
     .            string(i:i).eq.char(0) ) then
                  i = i - 1
              else
                  blanks = .false.
              end if
          else
*                                     ! Force exit at end of string
              blanks = .false.
          end if
      end do
 
***** Save length of string
      trimlen = i
 
***** Thats all
      return
      end
 
    

c********************************************************************

      integer function nblen (string)
c
c     given a character string, nblen returns the length of the string
c     to the last non-blank character, presuming the string is left-
c     justified, i.e. if string = '   xs  j   ', nblen = 8.
c
c     called non-library routines: none
c     language: standard fortran 77
c
      integer ls,i
      character*(*) string, blank*1, null*1
      data blank /' '/
c
      null = char(0)
      nblen = 0
      ls = len(string)
      if (ls .eq. 0) return
      do 1 i = ls, 1, -1
         if (string(i:i) .ne. blank .and. string(i:i) .ne. null) go to 2
    1    continue
      return
    2 nblen = i
      return
      end
c
c

c*********************************************************************            

      SUBROUTINE MONDAY(IDOY,IM,IDAY,IYR)
C
C     CONVERT A DAY OF YEAR TO MONTH AND DAY
C     RICK ABBOT - NOVEMBER 1984
C
      integer marray(12),idoy,iday,im,iyr,isum

      DATA MARRAY/31,28,31,30,31,30,31,31,30,31,30,31/
C
c     reset this to avoid problem with previous calls
      marray(2) = 28
      IF (MOD(IYR,4).EQ.0) MARRAY(2)=29
      ISUM=0
      DO 12 IM=1,12
      ISUM=ISUM+MARRAY(IM)
      IF (IDOY.GT.ISUM) GO TO 12
      IDAY=IDOY-ISUM+MARRAY(IM)
      GO TO 14
 12   CONTINUE
C
 14   RETURN
      END

c*************************************************************************

      Subroutine fix_y2k(iyear)

c     Check for 2 or 4-digit years.  If 2-digits, convert to 4 by guessing
c     the date.   Same as check_y2k except no warning printed.  R. King July 1999

      implicit none

      integer*4 iyear

      if( iyear.le.1900) then 
                  
c        earliest GPS launch date is 1978; earliest space-geodetic data about 1960 
         if( iyear.gt.60 ) then
            iyear = iyear + 1900
         else
            iyear = iyear + 2000
         endif  

      endif

      return
      end
      

