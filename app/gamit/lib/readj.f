      Subroutine READJ ( iuj,ischan,nchan,iprn,jd,t,svdt,icall,
     .                   jdtoc,toc,svcepc,svcrat,svcacc,valid)
C
C     Read SV clock corrections from an external data set (J-file)
C     R.W. King   21 December 1989 - from sb POLRED
c     K. Feigl    19 June, 1990      add svcepc,svcrat,svcacc, change name to READJ
c     R. King     20 July 1994 - change internal time from UTC to GPST

c     J-file specifications:
c
c     time tags are Satellite send times, valid on the interval centered on
c     that time tag, i.e. the record tagged at t(i) is valid from t1 to t2, where:
c                t1 = t(i) - (t(i)-t(i-1))/2
c                t2 = t(i) + (t(i+1)-t(i))/2
c   Input:
c     icall = 0    Initial call in makex/makex.f or model/setup.f to read in tabular values
c     icall = 1    Subsequent calls from makex.f or model.f to obtain clock values
c     jd           Julian day (GPST) at which to evaluate the polynomial.
c     t            Second of day (GPST) at which to evaluate the polynomial.
c     ischan       Array of satellite IDs (PRNs)
c     nchan        number of channels in array ischan
c     iprn         PRN number for which the SVDT is desired

c   Output 
c     jdtoc  Julian day (GPST) for the ref. time for polynomial expansion.
c     toc    Second of day (GPST) for the ref. time for polynomial expansion.
c     svdt   The value of the SV clock polynomial evaluated at (JD,T)
c     svcepc Sat. Vehicle. Clock Epoch offset in seconds
c     svcrat Sat. Vehicle. Clock Rate in dimensionless units
c     valid  This estimate should be good for 'valid' seconds after jdtoc,toc

      implicit none

      include '../includes/dimpar.h'  
      include '../includes/makex.h'
                                       
      character*1 gnss 
      character*3 j_time_flag
      character*16 fname
      character*60   afmt
      character*80  prog_name
      character*256  message,line

      integer*4 iuj,iprn,jd,ichan,icall,iyr,ihr,min,jd1
     .  ,julday,idoy,wkno,nchan,i,jdtoc,jprn,ioerr,
     .  jdsv(maxepc,maxsat),nrec(maxsat),ischan(maxsat),
     .  nblen,inqerr,len,rcpar

      real*8 t,svdt,xeaf0,xeaf1,xeaf2,tdiff,tdiffs,sec,sow,
     .       toc,svcepc,svcrat,svcacc,valid,sod,hrdiff,utcoff,taiutc,
     .       tsv(maxepc,maxsat),
     .       svcof1(maxepc,maxsat),
     .       svcof2(maxepc,maxsat),
     .       svcof3(maxepc,maxsat)


      save jdsv,tsv,svcof1,svcof2,svcof3,nrec

c     get the module and file names for report_stat calls
      len = rcpar(0,prog_name)  
      inquire ( unit=iuj,name=fname,iostat=inqerr )

c     maximum number of SV clock values for any one satellite

      if (icall .ne. 0 ) goto 200

C     Read the header records
                 
      read(iuj,'(1x)',iostat=ioerr) 
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/readj'
     .   ,fname,'Error reading first line of J-file',ioerr)
      read(iuj,'(a)',iostat=ioerr) line 
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/readj'
     .   ,fname,'Error reading second line of J-file',ioerr)
c     need the time-tag flag; can be in columns 18-20 or 21-23 
      if( line(18:20).eq.'UTC' .or. line(18:20).eq.'GPS' ) then
        read(line(18:20),'(a3)',iostat=ioerr) j_time_flag
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/readj'
     .   ,fname,'Error reading j_time_flag',ioerr)
      elseif ( line(21:23).eq.'UTC' .or. line(21:23).eq.'GPS' ) then
        read(line(18:20),'(a3)',iostat=ioerr) j_time_flag
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/readj'
     .   ,fname,'Error reading j_time_flag',ioerr)
      else
        call report_stat('FATAL',prog_name,'lib/readj'
     .   ,' ','Cannot find j_time_flag',ioerr)
      endif
      read(iuj,'(a60)',iostat=ioerr) afmt 
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/readj'
     .   ,fname,'Error reading format line of J-file',ioerr)
c     find the end of the format statement
      i = nblen(afmt)
c     choke if you can't find a format statement
      if (i .lt. 1 .or. afmt(1:1).ne.'(') 
     .   call report_stat('FATAL',prog_name,'lib/readj'
     .      ,fname,'Bogus format statement in J-file',ioerr)
      afmt = afmt(1:i)
      do i=1,maxsat
         nrec(i)= 0
      enddo

C     Read the data records into storage

100   continue   
c     use the presence of 'a1' in the format statement to determine if the
c     GNSS code is in the data record
      if( afmt(42:43).ne.'a1' ) then 
        read(iuj,fmt=afmt,end=120,iostat=ioerr)
     .     iyr,idoy,ihr,min,sec,wkno,sow,jprn,xeaf0,xeaf1,xeaf2   
      else
        read(iuj,fmt=afmt,end=120,iostat=ioerr)
     .     iyr,idoy,ihr,min,sec,wkno,sow,gnss,jprn,xeaf0,xeaf1,xeaf2   
      endif
      if( ioerr.ne.0 .and. ioerr.ne.-1 ) 
     .  call report_stat('FATAL',prog_name,'lib/readj'
     .              ,fname,'Error reading data record of J-file',ioerr)
      call fix_y2k(iyr)

      do 115 i=1,nchan
         if( jprn.eq.ischan(i) ) then
            if (nrec(i).ge.maxepc) then
               call report_stat('FATAL',prog_name,'lib/readj',' '
     .  ,'Number of clock coefficients in the J-file exceeds dimensions'
     .  ,0)
            else
               nrec(i) = nrec(i) + 1
            endif
c           Get julian day by adding idoy to JD(Jan. 1)
            jd1 = julday(1,1,iyr) + idoy - 1
c           second of day
            sod =  3600.d0*ihr+60.d0*min+sec
c           convert from UTC to GPST if necessary
            if( j_time_flag.ne.'GPS' ) then
              utcoff = taiutc(jd1) - 19.d0
              call timinc(jd1,sod,utcoff)
            endif
c           seconds after last record
            if (nrec(i) .gt. 1) then
               svdt  = 86400.d0*(jd1-jdsv(nrec(i)-1,i))
     .             + sod-tsv(nrec(i)-1,i)
            else
               svdt  = 10.d0
            endif

c           check to make sure that record is not duplicate:
            if (svdt .gt. 1.d0) then
               jdsv(nrec(i),i) = jd1
               tsv(nrec(i),i)  = sod

               svcof1(nrec(i),i) = xeaf0
               svcof2(nrec(i),i) = xeaf1
               svcof3(nrec(i),i) = xeaf2
            else
               nrec(i) = nrec(i) - 1
c**    Suppress this warning message - rwk 91/34/9
c               write (*,*) 'READJ: record for PRN '
c     .         ,ischan(i),' duplicates record ',nrec(i)
            endif
         endif
115   continue
      if (ioerr .eq. 0) goto 100

c     Check to make sure every satellite has a value
c     Also write a large time to the last time tag
120   continue
      do 125 i=1,nchan
         if( nrec(i).eq.0 ) then
           write(message,121) i,ischan(i)
121        format('Satellite channel ',i2,' (PRN=',i2,') not in J-file')
           call report_stat('FATAL',prog_name,'lib/readj',' ',message,0)
         endif
         tsv (nrec(i)+1,i) = tsv(nrec(i),i) + 1.0d0
         jdsv(nrec(i)+1,i) = jdsv(nrec(i),i) + 7
125   continue
      close (iuj)
      return

c------come here if ICALL>0--(not initial call)--------------------------------


200   continue
c     determine the channel
      ichan = 0
      do 204 i = 1,maxsat
         if (ischan(i) .eq. iprn) ichan = i
 204  continue
                  
      if (ichan .eq. 0) then
         write (message,'(a,i3,a)') 'READJ: PRN ',iprn
     .         ,' not found in J-file'
         call report_stat('FATAL',prog_name,'lib/readj',' ',message,0)
      endif

      if (nrec(ichan) .eq. 0) then
         write (message,'(a,i3)') 'No records in J-file for PRN ',iprn
         call report_stat('FATAL',prog_name,'lib/readj',' ',message,0)
      endif

C     Compute the clock correction using the closest value
      do 220 i = 1,nrec(ichan)
c        time of current record
         tdiff = (jd-jdsv(i,ichan))*86400.d0 + (t-tsv(i,ichan))
CD        write (*,*) jd,jdsv(i,ichan),t,tsv(i,ichan)
c        time of next record
         tdiffs = (jd-jdsv(i+1,ichan))*86400.d0 + (t-tsv(i+1,ichan))

c        check values until the difference starts increasing
         if( dabs(tdiff) .le. dabs(tdiffs)) then
            hrdiff = tdiff/3600.d0
            if( dabs(hrdiff).gt.24.04d0 ) then
               if( dabs(hrdiff).ge.100.d0 ) then
                 write(message,210) ischan(ichan),hrdiff
210              format('Clock coeff. for PRN '
     .                ,i2,' more than 100 hrs from obs. epoch ',
     .                ' (tdiff=',d12.3,' hrs)')
                 call report_stat('WARNING',prog_name,'lib/readj'
     .                           ,' ',message,0)
               else
                 write(message,211) ischan(ichan),hrdiff
211              format('Warning clock coeff. for PRN '
     .                ,i2,' more than 24 hrs from obs. epoch ',
     .                ' (tdiff=',f8.3,' hrs)')
                 call report_stat('WARNING',prog_name,'lib/readj'
     .                            ,' ',message,0)
               endif
            endif
            svdt =  svcof1(i,ichan)
     .          + svcof2(i,ichan)*tdiff
     .          + svcof3(i,ichan)*tdiff**2
            svcepc = svcof1(i,ichan)
            svcrat = svcof2(i,ichan)
            svcacc = svcof3(i,ichan)
            toc    = tsv (i,ichan)
            jdtoc  = jdsv(i,ichan)
c           estimate is valid halfway to the next one
            valid  = (86400.0d0*(jdsv(i+1,ichan)
     .                          -jdsv(i  ,ichan))
     .                           +tsv(i+1,ichan)
     .                           -tsv(i  ,ichan))/2.d0
            if (valid .lt. 1.d0) then
               write(message,'(a,f9.1,i3,2f9.1)') 'Invalid time tag: '
     .                                          , valid,jprn,t,toc
               call report_stat('FATAL',prog_name,'lib/readj',' '
     .                          ,message,0)
            else
               return
            endif
         endif
220   continue

c     come here if no entry
      write (message,'(a,i2,a,i2)')
     .  'SVDT value not selected for channel ',ichan
     .   ,' PRN ',ischan(ichan)
      call report_stat('WARNING',prog_name,'lib/readj',' ',message,0)
      valid  = 0.0d0
      svcepc = 0.0d0
      svcrat = 0.0d0
      svcacc = 0.0d0

      return
      end
