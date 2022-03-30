      Program ORBDIF
c
c Purpose:	Compare two orbits in t-file, SP1 or SP3 NGS format in terms 
c      of positional errors, either in x,y,z, radial, along, cross
c		track, or N,E,U components. Interpolation is made if the epoches
c		of the two orbits do not match. The tabuleted orbit with
c		longer time interval will automatically be selected as
c		reference orbit. If the time intervals are identical, then
c		the first is used as reference. Optionally a set of files 
c		containing delta_x/y/z or delta_radial/along/cross-track, or 
c      delta_vx/vy/vz, or delta_N/E/U and ground track position may
c		be generated for viewing/plotting purpose.
c
c P. Fang, pfang@pgga.ucsd.edu, April 1993
c S. McClusky simon@chandler.mit.edu September 1994
c
c input: Standard GAMIT T-Files
c		  NGS orbit A, in either sp1 or sp3 format
c		  NGS orbit B, in either sp1 or sp3 format
c output:	RMS summary of the overlap
c		a set of files (optional) that contain delta_x/y/z or
c		delta_radial/along/cross-track, or delta_vx/vy/vz, or 
c      delta_N/E/U and ground track position.
c
c variables: orbtyp   orbit type
c            ftype     file type T-file or NGS file 
c	          spfmt	    sp-file format (sp1 or sp3)
c		      org	    organization code (e.g. MIT, SIO)
c		      crdsys    coordinate system (e.g. ITR91)
c		      mjd	    integer part of starting M-julian day
c		      nprn	    number of satellites
c		      nrec	    number of epoches
c		      iprn	    satellite PRN numbers, set to 85 due to SP3
c				        format definition
c		      fmjd		fractional part of modified julian day
c		      stat/endt	stating/ending day in modified julian day
c		      delt		epoch time interval
c  orb(orbit_A/orbit_B interpolated to orbit A tab interval,
c  number of satellites,
c  number of epoches,
c	x/y/z/xdot,ydot,zdot components)
c
c note: maxepc has been set to the number of epoches of 7 day arc at
c		15 min tabular time interval		
c
      implicit none
c
      include '../includes/dimpar.h'
c
      CHARACTER* 1  ftype,frame,gnss
      character* 3  orbtyp(2),spfmt(2) 
      character* 4  icsnam(maxorb),org(2)
      character* 5  crdsys(2),precmod,nutmod,gravmod,i_frame,srpmod
     .           ,  eradmod,antradmod                        
      character*16  satnam(maxsat)
      character*24  tfname(2)
      character*120  version
      CHARACTER*80  orbnama,orbnamb,rmsnam 

      integer*4     luo(2),lunit(2),lurms,maxepc
      integer*4     mjd(2),nprn(2),nrec(2),iprn(2,85)
      integer*4     itsat(maxsat)
      integer*4     nics(2),iref,i2nd,iyrs,idoys,iyre,idoye
      integer*4     nintrs(2),nblen,namelen
      integer*4     iscrn,j,k,ksat1,ksat2,iepoch
      integer*4     ji01,jil1,jnow1,jlast1,iendf1,iy11,iy21
      integer*4     ji02,jil2,jnow2,jlast2,iendf2,iy12,iy22
      integer*4     iprnt,jdb,jdstp,jde
      integer*4     ierr  ! IOSTAT error on close (for temp file delete)

      parameter     (maxepc=673)
 
      real*8        fmjd(2),stat(2),endt(2),delt(2),statim,endtim
      real*8        orb(4,maxsat,maxepc,6),hrs,hre,hrt
      real*8        satics(maxorb,maxsat)
      real*8        deltma,deltmb,stat2ref,refepoch
      real*8        rvec1(maxytp),ytrp1(5,2,maxyt2,maxsat),rvec2(maxytp)
      real*8        yy1(10,maxytp,maxsat),statmjd,statfmjd,endmjd
      real*8        endfmjd,dnumint,truniref,truni2nd,timeiref,timei2nd
      real*8        yy2(10,maxytp,maxsat),ytrp2(5,2,maxyt2,maxsat)
      real*8        overlaps,overlape 
      real*8        tb,tstp,te

      logical  use_runstring  ! Set true if input from runstring

* MOD TAH 180106: Added reading from runstring
      

      data luo(1)/11/,luo(2)/12/,lurms/13/,iscrn/6/,iprnt/0/  
c
c     Get version number
c
      call oversn(version)
c
c     Get input output file names
* MOD TAH 201020: Added string initialization
      rmsnam = ' '
c

****  See if runstring passed
      call rcpar(1,ftype)
      if( ftype.ne.'T' .and. ftype.ne.'t' .and. 
     .    ftype.ne.'N' .and. ftype.ne.'n'  ) then   ! Use interactive
         use_runstring = .false.

         write(*,120)
 120     format('ORBDIF: Interactive input',/,
     .          'Runsting: <Type> <Frame> <T_file1> <T_File2> ',
     .          '<RMS File> <Plot> <Plot Ext>',/,
     .          'Differences are T_File2-T_file1')

         print *, ('Enter orbit file type T-Tile (T) or NGS (N) ')
         read(*, '(a)') ftype
* MOD TAH 200714: Allow exit with X
         if(ftype.eq.'X'.or.ftype.eq.'x') stop 'X option used'

         if (ftype.eq.'T'.or.ftype.eq.'t') then
           print *, ('Are T-Files Inertial (I) or Earth Fixed (E) ')
           read(*, '(a)') frame
           if (frame.eq.'I'.or.frame.eq.'i') then 
             i_frame = 'INERT'
           else if (frame.eq.'E'.or.frame.eq.'e') then 
             i_frame = 'EFIXD'
           else 
             write(iscrn,5) frame
5            format(1x,'incorrect reference frame specification ',a1)
             stop
           endif
           print *, ('Enter name of T-File (A) ')
	   read(*, '(a)') orbnama
	   open(unit=luo(1),file=orbnama,status='old',form='unformatted'
     .        , err=999)
	   print *, ('Enter name of T-File (B) ')
	   read(*, '(a)') orbnamb
	   open(unit=luo(2),file=orbnamb,status='old',form='unformatted'
     .        , err=999)
	   print *, ('Enter RMS output file name ')
	   read(*, '(a)') rmsnam
	   open(unit=lurms,file=rmsnam,status='unknown')
         else if (ftype.eq.'N'.or.ftype.eq.'n') then 
           write(iscrn,7)
7          format(53('*'),/
     .        ,'* Warning: When comparing orbits in an Earth fixed  *',/
     .        ,'* reference frame errors in calculation of radial   *',/
     .        ,'* along and cross track components will occur since *',/
     .        ,'* Earth fixed velocities are used in calculations   *',/
     .        , 53('*'),/)      
           print *, ('Enter NGS orbit file name A ')
	   read(*, '(a)') orbnama
	   open(unit=luo(1),file=orbnama,status='old')
	   print *, ('Enter NGS orbit file name B ')
	   read(*, '(a)') orbnamb
	   open(unit=luo(2),file=orbnamb,status='old')
	   print *, ('Enter RMS output file name ')
	   read(*, '(a)') rmsnam
	   open(unit=lurms,file=rmsnam,status='unknown')
         else
           write(iscrn,10)ftype
10         format(1x, 'Incorrect file type specification ',a1)
           stop       
         endif
      else    ! Read values from runstring.
         use_runstring = .true.
         call rcpar(2,frame)
         if (frame.eq.'I'.or.frame.eq.'i') then 
           i_frame = 'INERT'
         else if (frame.eq.'E'.or.frame.eq.'e') then 
           i_frame = 'EFIXD'
         else 
           write(iscrn,5) frame
           stop 'incorrect reference frame specification '
         endif
         call rcpar(3,orbnama)
         call rcpar(4,orbnamb)
* MOD TAH 200331: Open ASCII SP3 normally, unformatted 
*        open for t-files.  This seems a long standing
*        bug.  Not sure it ever worked for SP3 files.
         if( ftype.eq.'N'.or.ftype.eq.'n') then
	    open(unit=luo(1),file=orbnama,status='old', err=999)
	    open(unit=luo(2),file=orbnamb,status='old', err=999)
         else    ! t-file open unformatted (old code)
	    open(unit=luo(1),file=orbnama,status='old',
     .           form='unformatted' , err=999)
	    open(unit=luo(2),file=orbnamb,status='old',
     .           form='unformatted', err=999)
         endif
         call rcpar(5,rmsnam)
	 open(unit=lurms,file=rmsnam,status='unknown')
         if (ftype.eq.'N'.or.ftype.eq.'n') then 
           write(iscrn,7)
         endif 
      endif 


c
c read in the headers and get starting/ending times in M-julian day
c
        if (ftype.eq.'N'.or.ftype.eq.'n') then
	  call rdsphd(luo(1),mjd,fmjd,iprn,1,
     *          delt,nrec,nprn,orbtyp,org,crdsys,spfmt)
	  stat(1)=mjd(1)+2400001+fmjd(1)
	  endt(1)=stat(1)+((delt(1)*nrec(1))/86400.d0)
          if(orbtyp(1).eq.'e'.or.orbtyp(1).eq.'E') i_frame = 'EFIXD'
c          print*, 'stat(1) endt(1) ',stat(1),endt(1)
	
	  call rdsphd(luo(2),mjd,fmjd,iprn,2,
     *          delt,nrec,nprn,orbtyp,org,crdsys,spfmt)
          if(orbtyp(2).ne.orbtyp(1)) then
             print*, 
     .        'Warning NGS files are not in the same reference frame'
             print*, 'Frames are ',orbtyp(2), ' and ',orbtyp(1)
          endif 
	  stat(2)=mjd(2)+2400001+fmjd(2)
	  endt(2)=stat(2)+((delt(2)*nrec(2))/86400.d0)
c          print*, 'stat(2) endt(2) ',stat(2),endt(2)

        elseif (ftype.eq.'T'.or.ftype.eq.'t') then 
          call thdred( luo(1),iscrn,iprnt,nprn(1),gnss,itsat,satnam
     .               , jdb,tb,jdstp,tstp,delt(1),nrec(1),jde,te
     .               , nics(1),satics,nintrs(1),icsnam
     .               , precmod,nutmod,gravmod,i_frame,srpmod,eradmod
     .               , antradmod )

          do j = 1,nprn(1)
            iprn(1,j) = itsat(j)
          enddo   
          stat(1)=jdb + tb/86400.d0
          endt(1)=stat(1)+((delt(1)*nrec(1))/86400.d0)
          call thdred( luo(2),iscrn,iprnt,nprn(2),gnss,itsat,satnam
     .               , jdb,tb,jdstp,tstp,delt(2),nrec(2),jde,te
     .               , nics(2),satics,nintrs(2),icsnam
     .               , precmod,nutmod,gravmod,i_frame,srpmod,eradmod
     .               , antradmod )
          do j = 1,nprn(2)
            iprn(2,j) = itsat(j)
          enddo 
          stat(2)=jdb + tb/86400.d0
          endt(2)=stat(2)+((delt(2)*nrec(2))/86400.d0)
        else
          write(iscrn,100)ftype
100       format(1x, 'Incorrect file type specification ',a1)
          stop       
        endif
c          print*, 'stat(1) endt(1) ',stat(1),endt(1)
c          print*, 'stat(2) endt(2) ',stat(2),endt(2)
c
c compute the overlaped time range
c
	statim=max(stat(1),stat(2))
        statmjd=dint(statim)
        statfmjd=(statim-statmjd)*86400.d0
	endtim=min(endt(1),endt(2))
        endmjd=dint(endtim)
        endfmjd=(endtim-endmjd)*86400.d0   
c        print*, 'statim ',statim
c        print*, 'endtim ',endtim
c
c set reference orbit to one with shorter tabular interval T-file interpolation
c
          iref=1
          i2nd=2
	  if (delt(2).lt.delt(1)) then
		iref=2
		i2nd=1
          endif
c          print*, 'iref i2nd ',iref,i2nd
c
c read in both NGS orbits and store in tmp1.tmp and tmp2.tmp binary t-file 
c format for interpolation.
c
c        print*, 'nrec(iref) nrec(i2nd)',nrec(iref),nrec(i2nd)
c        print*, 'delt(iref) delt(i2nd) ',delt(iref),delt(i2nd)
        if (ftype.eq.'N'.or.ftype.eq.'n') then
          namelen = nblen(orbnama) 
* MOD TAH 200714: Removed debug print    
          ! print *,'orbnama, namelen: ',orbnama(1:namelen),namelen    
          tfname(iref)='tmp1.t'//orbnama(1:namelen)
          tfname(i2nd)='tmp2.t'//orbnama(1:namelen)
          call ttemp(luo(iref),tfname(iref),spfmt,nprn,nrec,iref
     .              ,lunit(iref),nintrs)
          rewind(lunit(iref))
          luo(iref)=lunit(iref)
          call ttemp(luo(i2nd),tfname(i2nd),spfmt,nprn,nrec,i2nd
     .              ,lunit(i2nd),nintrs)
          rewind(lunit(i2nd))
          luo(i2nd)=lunit(i2nd)
        endif
c        print*, 'luo(iref) luo(i2nd) ',luo(iref),luo(i2nd)
c
c Interpolate (i2nd) file with longest tabular interval so that interpolated 
c points are coincident with reference orbit tab intervals. Interpolate
c ref orbit too to obtain velocity values at tabular points. Store coincident
c tabular points in array orb for orbit overlapped periods (endt-stat).
c Cannot interpolate within 5 epochs of tfile end so step over these points.
c 
         write(iscrn,105)
105      format(/,1x,'Interpolating Orbit files')      
         stat2ref = statim-stat(iref)
c         print*,'stat2ref ',stat2ref
         dnumint = dint(((stat2ref*86400.d0)/delt(iref))+5.d0)
c         print*,'dnumint ',dnumint                      
         refepoch = stat(iref) + ((dnumint*delt(iref))/86400.d0)
c        print*, 'refepoch',refepoch
         truni2nd = refepoch - stat(i2nd)
         truni2nd = dnint(truni2nd*86400.d0) 
         truniref = refepoch - stat(iref) 
         truniref = dnint(truniref*86400.d0)
         do 106 j=1,1000   
         if(truni2nd.lt.5*delt(i2nd).or.truniref.lt.5*delt(iref)) then 
           truni2nd = truni2nd + (delt(iref)) 
           truniref = truniref + (delt(iref))
         else
           goto 107 
         endif 
106      continue
107      continue
c         print*, 'truniref truni2nd ',truniref,truni2nd
         do 115 iepoch = 1,10000 
            if(iepoch.gt.1) then
              truniref = truniref + (delt(iref))
              truni2nd = truni2nd + (delt(iref))
            endif
c         print*, 'truniref truni2nd ',truniref,truni2nd
            if (mod(iepoch,20).eq.0) print*, 'interval',iepoch
            timeiref = stat(iref) + (truniref/86400.d0)
            timei2nd = stat(i2nd) + (truni2nd/86400.d0)
c            print*,'timeiref timei2nd ',timeiref,timei2nd
c            print*, 'endt(iref) endt(i2nd) ',endt(iref),endt(i2nd)
            if (timeiref.ge.(endt(i2nd)-(5.d0*delt(i2nd)/86400.d0))
     .      .or.timeiref.ge.(endt(iref)-(5.d0*delt(iref)/86400.d0))) 
     .      goto 117
            do 110 ksat1=1,nprn(i2nd)
               call gsatel(2,truni2nd,luo(i2nd),ksat1,rvec1,ytrp1,yy1
     .                  , nprn(i2nd),delt(i2nd),nintrs(i2nd),ji01,jil1
     .                  ,iy11,iy21,jlast1,jnow1,iendf1,nrec(i2nd))
               do k=1,6
                  orb(i2nd,ksat1,iepoch,k)=rvec1(k)*1000.d0
               enddo     
Cd             write(*,109)(rvec1(k),k=1,6)
CD 109            format(1x,3(1x,f13.4),3(1x,f13.7))
110         continue
c
            do 112 ksat2=1,nprn(iref)
               call gsatel(2,truniref,luo(iref),ksat2,rvec2,ytrp2,yy2
     .                  , nprn(iref),delt(iref),nintrs(iref),ji02,jil2
     .                  ,iy12,iy22,jlast2,jnow2,iendf2,nrec(iref))
               do k=1,6
                  orb(iref,ksat2,iepoch,k)=rvec2(k)*1000.d0
               enddo
c              write(*,109)(rvec2(k),k=1,6)                    
112         continue  
115      continue
117      print*, 'interval',(iepoch-1)
c
c prepare and write out the summary file header
c
	 deltma=delt(1)/60
	 deltmb=delt(2)/60
         overlaps = refepoch
         overlape = (refepoch+(((iepoch-1)*delt(iref))/86400.d0))
	 call dayjul(int(overlaps),iyrs,idoys)
         call dayjul(int(overlape),iyre,idoye)
         hrs=(overlaps-dint(overlaps))*24
	 hre=(overlape-dint(overlape))*24
         hrt=(overlape-overlaps)*24
	 write(lurms,910) orbnama(1:nblen(orbnama)),deltma,
     *	 orbnamb(1:nblen(orbnamb)),deltmb,
     *   iyrs,idoys,hrs,iyre,idoye,hre,hrt
 910     format("Positional Errors (RMS in meters) between Two Orbits:"
     *         ,/,a," ( int:",f5.2,"min) vs. ",a," ( int:"
     *         ,f5.2,"min)",/,"Overlapped Time Range:",2(i5,i4,f6.2),
     *         " total:",f6.2," hours")

         write(lurms,920)
 920     format(/,84('-'),/,
     .        'PRN     Total    delta-X   delta-Y   delta-Z',
     .        '   d-Radial  d-Along   d-Cross      d-3D',/,
     .            84('-'))
c
c compare the two orbits
c   
* MOD TAH 200331: Passed refepoch (PEPJD) in to routine so that
*        MJD can be output into residual files.  
         call orbrms(lurms,iref,i2nd,iepoch,nprn,iprn,orb,delt,
     .               use_runstring, refepoch)
         print *,'Normal end of ORBDIF'

* MOD TAH 200714: Delete temporary files on end.
* MOD TAH 200810: Only delete files if SP3 files were input (t-files
*        created temporarily in this case.
         if (ftype.eq.'N'.or.ftype.eq.'n') then
            close(luo(iref),iostat=ierr, status='delete')
         else
            close(luo(iref),iostat=ierr)
         endif

         if( ierr.ne.0 ) write(*,930) ierr, iref, luo(iref)
 930     format('IOSTAT error ',i5,' Close/Deleting T-file ',i2,
     .          ' Unit ',i4)
* MOD TAH 200810: Only delete files if SP3 files were input (t-files
*        created temporarily in this case (2nd file; TAH 210201).
         if (ftype.eq.'N'.or.ftype.eq.'n') then
            close(luo(i2nd),iostat=ierr, status='delete')
         else
            close(luo(i2nd),iostat=ierr)
         endif

         if( ierr.ne.0 ) write(*,930) ierr, i2nd, luo(i2nd)
c
c suicide area
c
        goto 1100 
999     write(iscrn,1000)
1000    format(1x,'error opening t-file')
        stop 

1100	end
	subroutine rdsphd(lusp,mjd,fmjd,iprn,isw,
     *		delt,nrec,nprn,orbtyp,org,crdsys,spfmt)
c
c purpose:	read in sp1 or sp3 header information
c
	character*3 orbtyp(2),spfmt(2)
	character*4 org(2)
	character*5 crdsys(2)
	integer lusp,mjd(2),nrec(2),nprn(2),iprn(2,85),isw,iyr,mon
     .  , iday,ihr,min,igwk,i,ierr
	real*8 sec,delt(2),fmjd(2),dumy

c check the format
* MOD TAH 200714: Removed debug output
        ! print *,'Unit ', lusp, isw 
	read(lusp,'(a)',iostat=ierr) spfmt(isw)
        !print *, 'IOSTAT ',ioerr, spfmt(isw)
	read(lusp,'(a)',iostat=ierr) spfmt(isw)
        !print *, 'IOSTAT ',ioerr, spfmt(isw)
	rewind (lusp)

	if (spfmt(isw)(1:2).ne.'##') then
c
c read in sp1 format                                    
      org = ' ' 
      read(lusp,900) iyr,mon,iday,ihr,min,sec,delt(isw),mjd(isw),
     *		fmjd(isw),nrec(isw),orbtyp(isw),org(isw)(1:3)
900	format(4x,i4,4i3,f11.7,f15.7,i6,f16.14,1x,i6,a1,1x,a3)
	read(lusp,910) nprn(isw),(iprn(isw,i),i=1,34),
     *		crdsys(isw)(4:5),igwk
910	format(3x,i2,1x,34i2,a2,i4)
	spfmt(isw) = 'sp1'

	else

c read in sp3 format
	read(lusp,920) iyr,mon,iday,ihr,min,sec,nrec(isw),
     *		crdsys(isw),orbtyp(isw),org(isw)
920	format(3x,i4,4i3,1x,f11.8,2x,i6,7x,a5,1x,a3,1x,a4)
	read(lusp,930) igwk,dumy,delt(isw),mjd(isw),fmjd(isw)
930	format(3x,i4,1x,f15.8,1x,f14.8,1x,i5,1x,f15.13)
	read(lusp,940) nprn(isw),(iprn(isw,i),i=1,85)
* MOD TAH 060201: Changed 17i3 to 17(1x,I2) to allow for G01 etc.
940	format(4x,i2,3x,17(1x,i2),4(/,9x,17(1x,i2)),15(/))
	spfmt(isw) = 'sp3'

	endif

c convert delt from second to day
c	delt(isw)=delt(isw)/3600/24

        return
        end




