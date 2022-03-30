c  Generate a file for input to GMT pshistogram to generate a histogram of velocity residuals.
c     R. King 17 February 2004; last modified 18 October 2007

c The calling sequence is 
c
c    vel2hist [vel-file] [out-file] [N/E/U/H]  [max sigma] [max res] [max ratio] [sign/absv] 
c
c            [value/ratio/sigma] [scalesig] [min sigma] [lon1] [lat1] [lon2] [lat2]

c     Required:

c      [vel-file]  GLOBK prt/org or vel file (reads adjustments, correlations, and sigmas)
c      [out-file]  the name of the output file, containing a single column of values   
c      [N/E/U/H]   signifies to extract north, east, up, or horizontal velocity magnitude   
c      [max sigma] maximum sigma to be considered (for the quantity given in N/E/U/H below)
c      [max res ]  maximum residual to be considered
c      [max ratio] maximum ratio res/sigma to be considered        

c    Optional:

c      [sign/absv] signifies to retain the sign (for N,E,U) or take the absolute magnitude (default absv)
c      [value/ratio] signifies to use the residual value, the value/sigma, or the sigma itself  (default ratio)  
c      [scalesig ] scale the sigmas by this value   
c      [min sigma] minimum sigma to consider, allows (with max sigma) taking a band) (default 0)
c      [lon1 lat1 lon2 lat2] coordinates of a box within which to consider residuals (default all in file)

c The output is a single-column file of ratios to be plotted by pshistogram

      implicit none

      character*1 comp,col1
      character*5 valtype
      character*8 arg,site
      character*50 velfile,histfile
      character*256 line

      integer iuvel,iuout,iuhist,iuwrms,iunrms
     .       ,iclarg,iarg,nblen,num,ioerr

      real maxsig,maxres,maxrat,long,lat,evel,nvel,uvel
     .   , evadj,nvadj,evsig,nvsig,rho,uvadj,uvsig,hvres,hvsig,res,sig
     .   , rat,sum2e,sum2n,sum2h,sum2u,evnrms,nvnrms,hvnrms,uvnrms
     .   , sum2ew,sum2nw,sum2hw,sum2uw,evwrms,nvwrms,hvwrms,uvwrms 
     .   , histval,scalesig,minsig,lon1,lat1,lon2,lat2

      logical startvel,endvel,labs

c     Function
      logical inbox
                          
c     Assign the units
      iuvel = 1
      iuout = 2 
      iuhist = 3
      iuwrms = 4
      iunrms = 5              

c Print help if no arguments

      iarg = iclarg(1,arg)
      if( iarg.eq.0 ) then
        print *,
     .  'vel2hist needs at least 8 arguments--see utils/vel2hist.f'      
        stop
      endif
           
c Print some stuff to the screen

      print *,' '
      print *,' Starting vel2hist '
      print *,' ' 
      print *,'Note:  Default is normalized residual magnitude'
      print *,' ' 

      print *,'**This is the test version from 110430, maybe not tested'
      print *,'**If a problem revert to vel2hist.f.081022 in active' 

c Get the run-string arguments

      iarg = iclarg(1,velfile)
      open(unit=iuvel,file=velfile,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        print *,'Error opening velfile ',velfile,ioerr
        stop
      else
        print *,'Opened velfile ',velfile
      endif

      iarg = iclarg(2,histfile)
      open(unit=iuhist,file=histfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        print*, 'Error opening outfile ',histfile,ioerr
      else
        print *,'Opened output file ',histfile  
      endif
             
      iarg = iclarg(3,arg)   
      if( iarg.gt.0 ) then
        read(arg(1:1),'(a)') comp
      else
        print *,'Missing component in runstring'
        stop
      endif    
cd      print *,'read comp ',comp

      iarg = iclarg(4,arg)
      read(arg,'(f8.0)') maxsig 

      iarg = iclarg(5,arg)
      read(arg,'(f8.0)') maxres  
cd      print *,'read maxres ',maxres

      iarg = iclarg(6,arg)
      read(arg,'(f8.0)') maxrat    
cd      print *,'read maxmat ',maxrat
                        
      labs = .true.
      iarg = iclarg(7,arg)
      if( iarg.gt.0 ) then 
        if(arg(1:1).eq.'S' .or.arg(1:1).eq.'s' ) labs = .false.
cd         print *,'read labs ',labs
      endif                                      

      valtype = 'ratio'
      iarg = iclarg(8,arg)
      if( iarg.gt.0 ) then 
        valtype = arg(1:5)   
cd        print *,'read valtype ',valtype
      endif  
           
      scalesig = 1.0
      iarg = iclarg(9,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') scalesig
      endif                                                

      minsig = 0.
      iarg = iclarg(10,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') minsig
      endif

      lon1 = 0.
      lon2 = 360.
      lat1 = -90.
      lat2 = 90.
      iarg = iclarg(11,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') lon1
      endif
      iarg = iclarg(12,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') lat1
      endif
      iarg = iclarg(13,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') lon2
      endif
      iarg = iclarg(14,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') lat2
      endif

c Print the input to the screen

      print *,'Maximum sigma (mm)  residual (mm)  ratio: '
     .       ,maxsig,maxres,maxrat
      print *,'Component : ',comp
      print *,'Value type: ' ,valtype
      print *,'Abs value: ',labs  
      print *,'Minimum sigma (mm) : ',minsig
      print *,'Box corners (lon/lat) : ',lon1,lat1,' ',lon2,lat2
      print *,' '

      
c Open the output print file (record, not used for histogram) and write the header

      open(unit=iuout,file='vel2hist.out',status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        print*, 'Error opening print file','vel2hist.out',ioerr
      else
        print *,'Opened print file ','vel2hist.out'
      endif
      write(iuout,'(a)') 'Print file vel2hist.out '
      write(iuout,'(4a)') ' velfile: ',velfile(1:nblen(velfile))
     .   ,'  histfile: ',histfile(1:nblen(histfile))
      write(iuout,'(a,a1,a,a5)')' component: ',comp
     .   ,'  value type: ',valtype
      write(iuout,'(a,f4.1,a,2f4.1,a,f4.1,a,f4.1,/,a,4f8.2)') 
     .       ' scalesig=',scalesig,'  min/max sig=',minsig,maxsig
     .      ,'  maxres=',maxres,'  maxrat=',maxrat   
     .      ,'  Lon/Lat Box=',lon1,lat1,lon2,lat2
      write(iuout,'(a)') '* --> value not used '
                                                 
     
c Read the values from the vel file and write the output file
          
c     initialize the sums  
      sum2e = 0.
      sum2n = 0.  
      sum2h = 0.
      sum2u = 0.   
      sum2ew = 0.
      sum2nw = 0.  
      sum2hw = 0. 
      sum2uw = 0.
      num = 0
c     find the start of the velocity summary 
      startvel = .false.
      do while (.not.startvel ) 
        read(iuvel,'(a)',iostat=ioerr) line
cd        print *,'startvel line ',line
        if( ioerr.ne.0 ) then
           print *,'Error finding SUMMARY VELOCITY ',ioerr
           stop      
        endif
        if( line(2:12).eq.'SUMMARY VEL' ) startvel = .true.
      enddo 
c     skip lines until find a valid longitude  
      startvel = .false.
      do while (.not.startvel ) 
        read(iuvel,'(a)',iostat=ioerr) line    
cd        print *,'ioerr,LINE ',ioerr,line
        read(line,'(a1,f6.0)',iostat=ioerr) col1,long   
cd        print *,'ioerr col1 long ',ioerr,col1,long
        if( ioerr.eq.-1 ) then
          print *,'Unexpected EOF reading header '
          stop
        elseif( ioerr.ne.0 .or. abs(long).gt.360.) then     
cd          print *,'testing continue ioerr long ',ioerr,long
c         failure to read the floating pt number correctly
          continue
        else       
          startvel = .true.  
          backspace(iuvel)
        endif
      enddo  
      print *,'end of header ',line
c     read and write the values            
      endvel = .false. 
      do while (.not.endvel )     
        line = '  ' 
        read(iuvel,'(a)',iostat=ioerr) line  
cd          print *,'vel lines ',line   
        if( ioerr.eq.-1 ) then 
          print *,'EOF encountered in velocity file '
          endvel = .true.  
          print *,'Set endvel = true at eof'
        elseif( ioerr.ne.0 ) then
          print *,'Error reading velocity line ',ioerr  
          stop      
        elseif( line(1:1).ne.' ' ) then    
c         ignore non-blank lines (comments)
          continue
        elseif( line(6:6).ne.'.' ) then
c         end of valid velocity lines
          endvel = .true.         
          print *,'Set endvel = true at dot '  
        else         
c         determine whether old or new format
          if(nblen(line).gt.112 ) then
            read(line,'(2f11.5,1x,4f8.2,2f8.2,f7.3,2x,f8.2,2f8.2,1x,a8)'
     .        ,iostat=ioerr)      
     .        long,lat,evel,nvel,evadj,nvadj,evsig,nvsig,rho
     .       ,uvel,uvadj,uvsig,site   
cd            print *,' evadj nvadj evsig nvsig '
cd     .          ,evadj,nvadj,evsig,nvsig  
          else
            read(line,'(2f9.3,1x,4f8.2,2f8.2,f7.3,2x,f8.2,2f8.2,1x,a8)'
     .        ,iostat=ioerr)      
     .        long,lat,evel,nvel,evadj,nvadj,evsig,nvsig,rho
     .       ,uvel,uvadj,uvsig,site   
cd            print *,' evadj nvadj evsig nvsig '
cd     .          ,evadj,nvadj,evsig,nvsig  
          endif
          if( ioerr.ne.0 ) then
             print *,'Error decoding line values ',ioerr
             stop
          else      
c           see if selection by latitude and longitude
            if( inbox(long,lat,lon1,lat1,lon2,lat2) ) then    
c             trap vectors with 0 magnitude 
              if( evadj.eq.0 .and. nvadj.eq.0 ) then
                hvres = 0.
                hvsig = (evsig+nvsig)/2. 
              else
                call vmag(evadj,nvadj,evsig,nvsig,rho,hvres,hvsig)  
              endif
              if( comp.eq.'H' ) then
                res = hvres
                sig = hvsig*scalesig
                rat = hvres/(hvsig*scalesig)      
cd                print *,' site ',site
cd               print *,'evadj nvadj evsig nvsig, rho '
cd     .             ,evadj,nvadj,evsig,nvsig
cd                print *,' comp H res sig rat ',res,sig,rat
              elseif( comp.eq.'E') then
                res = evadj
                sig = evsig*scalesig
                rat = evadj/(evsig*scalesig)
              elseif( comp.eq.'N') then  
                res = nvadj
                sig = nvsig*scalesig
                rat = nvadj/(nvsig*scalesig)    
              elseif( comp.eq.'U') then
                res = uvadj
                sig = uvsig*scalesig
                rat = uvadj/(uvsig*scalesig)
              else
                print *,'Unknown component ',comp
                stop
              endif       
              if( labs) then
                res = abs(res)
                rat = abs(rat)
              endif   
              if( valtype.eq.'ratio' ) then
                histval = res/sig
              elseif ( valtype.eq.'value' ) then
                histval = res                   
              elseif ( valtype.eq.'sigma' ) then
                histval = sig
              else
                print *,'**Invalid value type ',valtype
              endif
              if( sig.ge.minsig.and.sig.le.maxsig.and.abs(res).le.maxres
c** I don't know why I put the _GPS qualifier here since EQ-affected
c** sites can/should be included.  If there are duplicate, linked sites,
c** this should be handled in another way.
c**     .           abs(rat).le.maxrat.and.site(5:8).eq.'_GPS') then  
     .            .and.abs(rat).le.maxrat )  then  
                write(iuhist,'(f6.2)') histval
                num = num + 1 
                sum2e = sum2e + (evadj/evsig)**2   
                sum2n = sum2n + (nvadj/nvsig)**2   
                sum2h = sum2h + (hvres/hvsig)**2  
                sum2u = sum2u + (uvadj/uvsig)**2  
                sum2ew = sum2ew + 1./evsig**2
                sum2nw = sum2nw + 1./nvsig**2     
                sum2hw = sum2hw + 1./hvsig**2   
                sum2uw = sum2uw + 1./uvsig**2    
cd                print *,'num ev/es sum2e ',num,evadj/evsig,sum2e
                write(iuout,'(1x,a4,3f6.2)') site,res,sig,histval
              else
                write(iuout,'(a,a4,3f6.2)')  '*',site,res,sig,histval
c             endif on sigma test
              endif
c           endif on lat/lon box test
            endif
c         endif on valid decode of line         
          endif  
c       endif on reading line
        endif
      enddo 

c* Write the statistics to the screen, the output print file,
c* and short files 'wrms' and 'nrms' files to be read by sh_velhist

      evnrms = sqrt( sum2e/float(num) ) 
      nvnrms = sqrt( sum2n/float(num) )  
      hvnrms = sqrt( sum2h/float(num) ) 
      uvnrms = sqrt( sum2u/float(num) )  
      evwrms = sqrt(sum2e/sum2ew)
      nvwrms = sqrt(sum2n/sum2nw)   
      hvwrms = sqrt(sum2h/sum2hw)
      uvwrms = sqrt(sum2u/sum2uw)
      print *,num,' sites included '  
      write(iuout,'(1x,i5,a)') num,' sites included ' 
      print *,'WMRS E N H U: ',evwrms,nvwrms,hvwrms,uvwrms
      write(iuout,'(a,4f6.2)')  
     .  'WMRS E N H U: ',evwrms,nvwrms,hvwrms,uvwrms
       print *,'NRMS E N H U: ',evnrms,nvnrms,hvnrms,uvnrms 
      write(iuout,'(a,4f6.2)')     
     .  'NRMS E N H U: ',evnrms,nvnrms,hvnrms,uvnrms 
      print *,'Sigmas scaled by ',scalesig
     .    ,' for histogram but not NRMS'  
      write(iuout,'(a,f4.2,a)') 'Sigmas scaled by ',scalesig
     .    ,' for histogram but not NRMS' 
      open(unit=iuwrms,file='vel2hist.wrms',status='unknown'
     .    ,iostat=ioerr)
      if( ioerr.ne.0) print *,'Error opening vel2hist.wrms' 
      write(iuwrms,'(a,i4,a,4f6.2)') "N=",num,"   WRMS E N Hmag U "
     .     ,evwrms,nvwrms,hvwrms,uvwrms                
      open(unit=iunrms,file='vel2hist.nrms',status='unknown'
     .    ,iostat=ioerr)
      if( ioerr.ne.0) print *,'Error opening vel2hist.nrms' 
      write(iunrms,'(a,i4,a,4f6.2)') "N=",num,"   NRMS E N Hmag U  "
     .     ,evnrms,nvnrms,hvnrms,uvnrms       

      stop
      end

c-----------------------------------------------------------

      Subroutine vmag(e,n,esig,nsig,rho,v,vsig) 

      implicit none

      real e,n,esig,nsig,rho,v,vsig

      v = sqrt(e**2 + n**2)
      vsig = sqrt( (e*esig)**2 + (n*nsig)**2 +
     .                     2.*rho*e*esig*n*nsig ) / v  
      return
      end                                                    

c------------------------------------------------------------

      Function inbox( long,lat,lon1,lat1,lon2,lat2 ) 
                    
      implicit none
      logical inbox 
      real long,lat,lon1,lat1,lon2,lat2

      inbox = .true.
      if( long.lt.lon1 .or. long.gt.lon2  .or.
     .    lat.lt.lat1 .or. lat.gt.lat2 ) inbox = .false.     
      return
      end
      

c------------------------------------------------------------

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

c----------------------------------------------------------

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
 
c-------------------------------------------------------------------

CTITLE 'TRIMLEN'
 
 
      integer*4 function trimlen(string)
 
 
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
              if( string(i:i).eq.' ' ) then
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
 
        


                              


