Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.   All rights reserved.
C
      Subroutine arcinp(tname,gname,jds,ts,jdf,tf,stepsize)

c     Subroutine to create an ARC input file
c
c     INPUT:  tname  bctot created tfile name   C*16
c             gname  bctot created gfile name   C*16    
c             jds    session start PEP JD       I*4
c             ts     session start sec-of-day   R*8
c             jdf    session finish PEP JD      I*4
c             tf     session finish  sec-of-day R*8

c     S. McClusky 1995/10/11; R.King 2013/12/17 for new input-file format 
c     R. King 2016/8/4 change the time arguments
c
      implicit none
c
      include '../includes/dimpar.h'   

c
      integer*4 nics,nsat,nepchs,itsat(maxsat)
     .        , idoys,idoyf
     .        , iut,iscrn,iprnt,i,j,iac
     .        , jds,iyrs,hrs,mins,jdf,iyrf,hrf,minf
     .        , jdb,jde,jdstp,nintrs
     .        , ioerr

      real*8    ts,tf, doys,doyf,satics(maxorb,maxsat),secs,secf
     .        , tb,te,tstp,sdelt,stepsize

      character*1  gnss(maxsat)
      character*4  end,icsnam(maxorb),g_time_flag
      character*5  precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*10 arcin_file
      character*16 tname,gname,satnam(maxsat),buf16
      character*20 frame_name

c
      data end/'END '/
      data iut/16/,iac/18/,iscrn/6/,iprnt/0/
      arcin_file = 'bc_arc.inp'  
c
c        Open the T-file
      call lowers(tname)
      open(unit=iut,file=tname,status='old', form='unformatted'
     .,iostat=ioerr)
      if (ioerr .ne. 0 ) then
         call report_stat('FATAL','BCTOT','orbits/arcinp',tname,
     .    'Error opening T-file: ',ioerr)
      else
        call report_stat('STATUS','BCTOT','orbits/arcinp',tname,
     .    'Opened broadcast T-file: ',ioerr)
      endif
c
c        Open the arc input file
      open(unit=iac,file=arcin_file,status='unknown',form='formatted',
     .iostat=ioerr)
      if (ioerr .ne. 0 ) then
        call report_stat('FATAL','BCTOT','orbits/arcinp',arcin_file,
     .    'Error opening ARC input file: ',ioerr)
      else
        call report_stat('STTUS','BCTOT','orbits/arcinp',arcin_file,
     .    'Opened  ARC input file: ',ioerr)
      endif
c
c        Read the T-file header
      frame = 'J2000'
      call thdred ( iut,iscrn,iprnt,nsat,gnss,itsat,satnam
     .            , jdb,tb,jdstp,tstp,sdelt,nepchs,jde,te
     .            , nics,satics,nintrs,icsnam
     .            , precmod,nutmod,gravmod,frame,srpmod
     .            , eradmod,antradmod )
c
c        Write the ARC input batchfile
      call dayjul(jds,iyrs,idoys)
      call ds2hms(iyrs,idoys,ts,hrs,mins,secs)
      call dayjul(jdf,iyrf,idoyf)
      call ds2hms(iyrf,idoyf,tf,hrf,minf,secf)
      g_time_flag = 'GPST'
      do j=1,nsat
        write(iac,'(a16)') satnam(j)
      enddo
      write(iac,'(a4)') end  
      write(frame_name,'(a,a5)') 'INERTIAL       '//frame               
      write(iac,75) gravmod,srpmod,sdelt,stepsize,g_time_flag,frame_name
     .              ,precmod,eradmod,antradmod 
   75 format(a5,1x,a5,1x,f6.1,1x,f7.2,3x,a4,1x,a20,3(1x,a5))
      write(iac,80)
   80 format('arcout.bctot')
      write(iac,85)gname
   85 format(a10)
      write(iac,90)
   90 format(1x)
      write(iac,95) iyrs,idoys,hrs,mins,secs     
      write(iac,95) iyrf,idoyf,hrf,minf,secf
   95 format(i4,1x,i3,1x,i2,1x,i2,1x,f8.5)
      write(iac,105)
  105 format('Y')
      write(iac,110)
  110 format('tfile.tmp')
      close (iac)
          call report_stat('Status','BCTOT','orbits/arcinp',arcin_file,
     .    'Successfully wrote ARC-input file: ',0)
      return
      end
