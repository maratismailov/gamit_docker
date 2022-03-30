      subroutine sdtrit( iut,iungs,iuclkfile,iprnt,spfile,nintrs,nsat
     .                 , gnss,itsat,spfmt,otlmod
     .                 , jds,ts,sdelt,nepcht,nstart,nstop, sample_rate )
c
c       Read a T-file and write the data records in an NGS SP#1 or SP#3 file
c       R.W. King  November 1991
c	P. Fang modified to output SP$3 records, March 1993 
c  R. King modified to output SP3-C  September2009 
c  T.Herring modified to allow sample_rate to be passed. 210122
c
      implicit none
c
      include '../includes/dimpar.h'

      integer*4 iut,iungs,iuclkfile,iprnt,jd,jds,itsat,nepchs
     .        , nepcht,icount,isat,nsat,nintrs,i,j
     .        , iyr,idoy,imon,iday,ihr,imin
     .        , ji0,jil,iy1,iy2,jnow,jlast,iendf
     .        , nstart,nstop,iepcht,iepchs,notl

* MOD TAH 210122: Introduced higher sampling rate factor (e.g. for 900-sec sampled
*     t-file to output at 300-seconds, sample_rate = 3
      integer*4 sample_rate   ! Sampling rate multiplier on t-file rate;
                              ! Value must >= 1; default is 1.

      real*8 x,ts,sec,sdelt,sod,rvec,ytrp,yy,trun,svclk(maxsat),docmc(3)
     
      character*1 gnss
      character*2 buf2,ksat
      character*3 spfmt
      character*4 asv
      character*8 otlmod
      character*16 spfile
      character*256 message

      integer*4 maxtrm
      parameter( maxtrm=(3+maxorb)*maxsat)

      dimension itsat(maxsat),asv(maxsat),x(maxtrm,maxsat)
C
      DIMENSION RVEC(MAXYT2)
      DIMENSION YTRP(5,2,MAXYT2,MAXSAT),YY(10,MAXYTP,MAXSAT)
C
c Initialize the counters for T-file and NGS file epochs (checked at end)

      iepcht= 0   ! TAH 210122: This variable does not seem to do anything.
      iepchs= 0   

c Initialize the interpolation pointers

      JI0= 0
      JIL= 0
      IY1= 0
      IY2= 0
      JLAST= 0
      IENDF= 0

C Set the SV ids  (modified rwk 090925 for SP3-C)
     
      do i=1,nsat     
        asv(i) = '    '
        write(asv(i),'(2x,i2)') itsat(i)
        asv(i)(1:1) = 'P'    
        asv(i)(2:2) = gnss
        if( asv(i)(3:3).eq.' ') asv(i)(3:3) = '0'
      enddo


c Write the data records
       
       nepchs = nstop - nstart + 1  

* MOD TAH 210122: Adjust the start time          
C      trun = (nstart-2) * sdelt 
       trun = (nstart-1) * sdelt -  sdelt/sample_rate 
 
* MOD TAH 210122: Increased icount end by sample_rate factor.
       do 50 icount=1,nepchs*sample_rate

* MOD TAH 210122: Decreas the time-step by sample_rate factor
         trun = trun + sdelt/sample_rate
c        compute the time of the epoch
         jd= jds
         sod= ts
         call timinc( jd,sod,trun )
         call dayjul( jd,iyr,idoy )
         call monday( idoy,imon,iday,iyr )
c        avoid roundoff by assuming that the epochs are even seconds
         sod = anint(sod)
         call ds2hms( iyr,idoy,sod,ihr,imin,sec )     
c        get the CMC corrections for ocean tidal loading
         if( otlmod(1:4).ne.'NONE')  then              
c          hard wire # of components
           notl = 11
           call otlcmc( jd,sod,otlmod,notl,2,docmc )   
c*          write(6,'(a,i8,f8.1,3f7.4)') 'CMC db ',jd,sod,(docmc(i),i=1,3)
         endif
c        read or interpolate the SV coordinates and clock values
         do  isat = 1, nsat
           call gsatel( 2,trun,iut,isat,rvec,ytrp,yy,nsat,sdelt
     .                , nintrs,ji0,jil,iy1,iy2,jlast,jnow,iendf
     .                , nepcht )
c          check for unexpected end-of-file
           if( iendf.ne.0 ) then
              if( icount.eq.1 ) then
                 goto 110
              else
                 goto 120
              endif
           endif
           do  i=1,6
             x(i,isat) = rvec(i)
           enddo   
c          apply the CMC corrections for the ocean loading model  
           if( otlmod(1:4).ne.'NONE' ) then
             do i=1,3
               x(i,isat) = x(i,isat)- docmc(i)/1.d3
             enddo
           endif
         enddo     

c        get the clock values              
         if( iuclkfile.gt.0 ) then   
* MOD TAH 201211: Added gnss to calling arguments to get correct
*          gnss clock.  Routine only seemed to be called here.           
           call rdsvclk (iuclkfile,iprnt,jd,sod,nsat,gnss,itsat,svclk )
         else
           do i=1,nsat  
             svclk(i) = .999999999999d0
           enddo
         endif
  

         iepcht = iepcht + 1
                           
c        Write the SP3-C records

	if (spfmt.eq.'sp3'.or.spfmt.eq.'SP3') then

         write(iungs,30) iyr,imon,iday,ihr,imin,sec
   30    format('*  ',i4,4(1x,i2),1x,f11.8)
         iepchs = iepchs + 1
         do j=1,nsat
           write(iungs,35) asv(j),(x(i,j),i=1,3),svclk(j)*1.d6
   35      format(a4,3(1x,f13.6),1x,f13.6)
         enddo

	else
c        SP1 
         call report_stat('FATAL','TTONGS','orbits/sdtrit',' '
     .                   ,'Only SP3 format now supported',0)

	endif

   50  continue

* MOD TAH 201211: Add EOF to end of file
      write(iungs,'(a)') 'EOF'

  
      write(message,'(a,i5,a,a3,a,a16 )') 
     .        'Wrote ',iepchs,' epochs to ',spfmt,' file ',spfile           
      call report_stat('STATUS','TTONGS','orbits/sdtrit',' ',message,0)
      return

  110 write(message,111) iyr,imon,iday,ihr,imin,sec
  111 format('End of file while looking for starting time on T-file.'
     2      ,' Last epoch read: ',i4,4(1x,i2),1x,f10.7 )
      call report_stat('FATAL','TTONGS','orbits/sdtrit',' ',message,0)

  120 write(message,121) jds,ts,jd,sod
  121 format('Unexpected end of file in data records,'
     2      ,' Last computed: JDS,TS,JD,T ',2(i7,f8.5) )
      call report_stat('FATAL','TTONGS','orbits/sdtrit',' ',message,0)

      end
