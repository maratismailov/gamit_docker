      Program GTOG

c Purpose:  To take a gfile in whatever inertial frame and precession model
c           and rotate it to whatever frame and precession model is
c           requested. This is a program to overcome problems which may be
c           encountered when SIO routine processing moves to J2000 inertial
c           frame.
c
c NOTE:     Program uses subroutine ROT_GFILE from ARC to rotate the gfile
c
c Paul Tregoning 25th July, 1995   
c R. King 22 May 2012: replace call to lib/ghdred with one to /lib/read_gfile

      implicit none

      include '../includes/dimpar.h'

      integer i,iugin,iugout,ioerr,ngsats,nics,ngics,jde,jdb,jdf

      real*8 trun0,trunf,te,tb,tf,delt,tdtoff,satics(maxorb,maxsat)

      character*1 upperc
      character*4 time_type,icsnam(maxorb)
      character*5 f_in,f_out,pr_in,pr_out,ans,srpmod,nutmod,gravmod
     .          , eradmod,antradmod
c rwk      character*10 ingfile  
      character*16 satnam(maxsat),ingfile,outgfile
c rwk      character*22 outgfile
      character*120 version
      character*256 message

      logical moveold/.false./

c units - ut1. pole. nutabl.
       integer*4 iscrn,iterm,inut,iut1,ipole

      parameter (iugin = 10,iugout = 11,iterm = 5)
      ioerr = 0
      iscrn = 6

c write out a status line
      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,5)version
    5 format('Started GTOG ',a120)
      call report_stat('STATUS','GTOG','orbits/gtog',' ',message,0)

c  get the name of the inputand output gfiles
      ingfile = ' ' 
      outgfile = ' ' 
      write(iscrn,'(a)')' Enter input gfile name : '
      read(iterm,'(a)')ingfile
      write(iscrn,'(a)')' Enter output gfile name : '
      read(iterm,'(a)')outgfile

c open the eop files
      call open_eop(inut,iut1,ipole)

c  open the input and output gfiles
      open(iugin,file=ingfile,status='old',iostat=ioerr)
      if (ioerr .ne. 0 ) then
         call report_stat('FATAL','GTOG','orbits/gtog',ingfile,
     .   'Error opening input G-file: ',ioerr)
      endif
      open(iugout,file=outgfile,status='new',iostat=ioerr)
      if(ioerr.ne.0)then
        write(iscrn,10)outgfile
10      format(/,5x,'File ',a22,' exists. Overwrite (y/n) ? ')
        read(iterm,'(a1)')ans(1:1)

        if(ans(1:1).eq.'y'.or.ans(1:1).eq.'Y')then
          open(iugout,file=outgfile,status='unknown',iostat=ioerr)
          if (ioerr .ne. 0 ) then
            call report_stat('FATAL','GTOG','orbits/gtog',outgfile,
     .      'Error opening output G-file: ',ioerr)
          endif
        else
          call report_stat('STATUS','GTOG','orbits/gtog',outgfile,
     .    'Stopped in GTOG output G-file exists: ',ioerr)
        endif
      endif

c read the gfile header to get the gfile frame/precmod and common info
c      nics = 9
c      call ghdred ( ingfile,nics,f_in,pr_in,ngics,iugin,'GTOG'
c     .             ,trun0,trunf,te,tb,tf,jde,jdb,jdf,delt,tdtoff
c     .             ,time_type,srpmod,nutmod,gravmod,icsnam,satnam )

c read the gfile 
      call read_gfile( ingfile,iugin,jde,te,f_in,pr_in
     .               , nutmod,gravmod,srpmod,eradmod,antradmod
     .               , ngsats,nics,icsnam,satnam,satics ) 

      write(iscrn,'(a,a5,a,a5,/)')' Gfile Inertial frame - ',f_in
     .                          ,' Precession - ',pr_in

c  find out the final required inertial frame
20    write(iscrn,'(a)')' Enter required inertial frame (B1950/J2000):'
      read(iterm,'(a5)')f_out
      f_out(1:1) = upperc(f_out(1:1))

c set default precession model
      if(f_out.eq.'J2000')then
        pr_out = 'IAU76'
      elseif(f_out.eq.'B1950')then
        pr_out = 'IAU68'
      else
        write(iscrn,'(a,a)')' Invalid frame: ',f_out
        goto 20
      endif

c get precession model - show the default option
      write(iscrn,'(a,a5,a)')' Enter required precession model ('
     .                        ,pr_out,'): '
      read(iterm,'(a5)')ans

      if(ans.ne.' ')then
        do i=1,5
          pr_out(i:i) = upperc(ans(i:i))
        enddo
      endif

c  now call rot_gfile and rotate to the rquired inertial space
      call rot_gfile(ingfile,f_out,f_in,pr_out,pr_in,moveold,iugout
     .               ,iugin,'GTOG',te,jde
     .               ,tdtoff,inut,iut1,ipole )

      call report_stat('STATUS','GTOG','orbits/gtog',' ',
     .'Normal end to GTOG ',0)

      stop
      end
