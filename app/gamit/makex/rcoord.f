      subroutine rcoord  (lu,site,dlat,mlat,seclat,latflag,
     .   dlon,mlon,seclon,lonflag,radius )

c
c     Read a priori estimates of the site coordinates from
c     the coordinates file.
c
c     These coordinates are used to compute clock offsets
c     and to label the Xfile.
c
c     Generalize L-file to have N/S and E/W longitudes - Yehuda Bock 4/25/90
c
C     New conventions for coordinates are:
c        Right-handed coordinate system.
c        latitude can be  'N' or 'S'
C        longitude can be 'E' or 'W'
C        Yehuda Bock 4/25/90
c     For history's sake, the OLD convention was:
c        Left-handed coordinate system.
c        latitude could be  ' ' (positive) west  or '-' (negative) east
c        longitude could be ' ' (positive) north or '-' (negative) south
c        This routine is back compatible with this convention.
c
c     DATA DECLARATIONS
c
C     logical unit (assumed open)

      integer*4     lu
      integer*4     dlat,mlat,dlon,mlon,ios

      real*8        seclat,seclon,radius

      character*1 latflag,lonflag                                              
      character*4 site,site4,jscr4
      character*64 jscr
c     case conversion functions:
      character*1 lowerc


      include '../includes/makex.h'

      rewind(lu)

c     read a string from the coordinates file (L-file format)
  601 continue
      read (lu,602,end=604,err=604,iostat=ios) jscr
  602 format (a61)

c     compare string to site   
      jscr4 = jscr(1:4)
      call lowers(jscr4)
      site4 = site      
      call lowers(site4)
      if (jscr4.eq.site4 ) then
         read (jscr,3255)
     .   latflag,dlat,mlat,seclat,lonflag,dlon,mlon,seclon,radius
 3255    format(17x,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f12.3)

C        To accommodate earlier GAMIT assumption that positive
c        values in L-file indicate north latitude and west longitude
         if(latflag.eq.' ') latflag='N'
         if(lonflag.eq.' ') lonflag='W'
         if(latflag.eq.'-' .or. dlat.lt.0) latflag='S'
         if(lonflag.eq.'-' .or. dlon.lt.0) lonflag='W'
         latflag = lowerc(latflag)
         lonflag = lowerc(lonflag)
      else
c        Read another line
         goto 601
      endif
C
      return

C     come here on error or end of file
604   continue
      write (uscren,1604) site
      write (uinfor,1604) site
1604  format (//,' RCOORD: Station ',a4,
     .     ' not found or not readable in coordinates file.')
      call ferror(ios,uscren)
      call ferror(ios,uinfor)
C
      end
