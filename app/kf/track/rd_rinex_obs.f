CTITLE READ_RINEX_OBS
 
      subroutine read_rinex_obs_1(eof, buffer,in_code, curr_stopgo)

      implicit none
 
*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
*     Split this to 2 parts, 
*            read_rinex_obs_1 : read in date information
*            read_rinex_obs_2 : read in obs. data
*     Gang Chen  971011
      include '../includes/const_param.h'
      include '../includes/xfile_def.h'
      include 'track_com.h'

*   in_code     -  0 normal read from file; 
*               -  1  from self buffer
      integer*4 in_code
*   date(5)     - Ymdhm of observation
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5),  ierr, jerr, i,id, trimlen, len
      integer*4 nline, lent, k   ! used for reading multiple lines.
c      character*256 line, buffer
      character*80 line, buffer
*   vals(5)     - Range values read
      character*1 xf_svtype(xf_maxsat)    ! Type of satellite 
                    ! R=GLONASS; G=GPS, E=GALELEO (converted to blank)
 
      real*8 sectag
 
*   eof     - Indicates end of file.
*   curr_stopgo -- status: False when receiver static
*   epfound -- Set true when EPFOUND
 
      logical eof, curr_stopgo, epfound
      logical eoh  ! End of header, used when RX files are cat'd
      integer*4 indx   ! Position in string
      character*1 cr   ! Carriages return

      cr = char(13)
  
****  Read in the next line from file
      if(in_code.eq.0) then
          read(obs_lu,'(a80)', iostat=ierr) line
*         Clean the "G" in firset line   
* MOD TAH 041216: Not needed anymore        
!         call clean_char(line,"G", 256)
          write(buffer,'(a80)') line
      else

****      Read in the next line from buffer
          read(buffer,'(a80)') line  
          in_code = 0        
      endif
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if

*     Decode the line
      epfound = .false.
      do while ( .not.epfound )
         epfound = .true.  
         len=trimlen(line)

* MOD TAH 100501: Feature from trackRTr:  Allow concatenated RX files
*        Check to see if rinex header line (from concatinated
*        file).  If so then skip new header
         call casefold(line)
         indx = index(line,'OBSERVATION DATA')
         if( indx.gt.0 ) then
             write(*,'(a,1x,a)') 'New RINEX header found in',
     .             trim(obs_file(obs_lu-50))
*            Skip to end of header
             eoh = .false.
             do while ( .not.eoh .and. ierr.eq.0 )
                 read(obs_lu,'(a)', iostat=ierr) line 
                 call sub_char(line,cr,' ')
                 call casefold(line)
                 call report_error('IOSTAT',ierr,'read', 
     .                 obs_file(obs_lu-50),1,'End of header not found')
                 if( trimlen(line).lt.10) eoh = .true.
                 if( index(line,'END OF HEAD').gt.0 ) eoh = .true.
             end do
*            Now read next line
             read(obs_lu,'(a80)', iostat=ierr) line
         end if

*        Now decode line.
         read(line,210,iostat=ierr) date, sectag, id, xf_msat,
     .               (xf_svtype(i),xf_iprn(i),i=1,min(xf_msat,12))
 210     format(5i3,f11.7,i3,i3,24(a1,i2))
* MOD TAH 061220: Only read the extended line if the ID is not
*     between 2 and 4 (comment and stop/go records)
         if( (id.lt.2 .or. id.gt.4)) then
             nline =  int((xf_msat-1)/12)
             do k = 1, nline
                lent = 12
                if( k.eq.nline ) lent = xf_msat-nline*12 

                read(obs_lu,'(a)',iostat=ierr) line
                read(line,215,iostat=jerr) 
     .                    (xf_svtype(i+k*12),xf_iprn(i+k*12),i=1,lent)
                call report_error('IOSTAT',jerr,'decod',line,0,
     .                    'GNSS Type and PRN')
 215            format(32x,12(a1,i2))
             end do  
         end if    

         call ymdhms_to_jd( date, sectag, xf_jd)
         call ymdhms_to_mjd( date, sectag, xf_mjd)           
*        Check to see if power-fail and comment indicator
* MOD TAH 040316: Added checks on other modes: IDs are
*        2 - Start kinematic mode
*        3 - Start static mode
*        4 - Power fail 
         if( id.ge.2 .and. id.le.4 ) then
*            OK, Skip over the next group of lines
             do i = 1, xf_msat
                read(obs_lu,'(a)',iostat=ierr) line
                write(*,220) date,sectag,i,line(1:40)
 220            format('RX Comment at ',5i3,1x,f8.3,' # ',i3,1x,a)
             end do
* MOD TAH 0401316: Set curr_stopgo status is ID changes it
             if( id.eq.2 ) curr_stopgo = .true.
             if( id.eq.3 ) curr_stopgo = .false.

*            Now get the next observation record
             read(obs_lu,'(a80)', iostat=ierr) line
* MOD TAH 160412: Not needed.
!            call clean_char(line,"G", 256)
             epfound = .false.
             if( ierr.ne.0 ) epfound = .true.
          end if
      end do

*         Clean the "G" in firset line          
           
      if( ierr.ne.0 ) then
           eof = .true.
           if( ierr.ne.-1 ) then 
              call report_error('IOSTAT',ierr,'read',line,0,
     .                          'RX_RINEX_OBS_1')
              write(*,*) " Warning rx_rinex_obs_1 error msg ", ierr
              write(*,*) " Make sure Rinex is unix file ", ierr
           end if
      endif

****  Convert PRN number of non-GPS satellites to values greater than 32
      do i = 1, xf_msat
         if( xf_svtype(i).eq.'R' ) then
             xf_iprn(i) = xf_iprn(i) + 100
         end if
         if( xf_svtype(i).eq.'E' ) then
             xf_iprn(i) = xf_iprn(i) + 200
         end if
         if( xf_svtype(i).eq.'C' ) then
             xf_iprn(i) = xf_iprn(i)+  300
         end if
         if( xf_svtype(i).eq.'J' ) then
             xf_iprn(i) = xf_iprn(i) + 400
         end if
         if( xf_svtype(i).eq.'S' ) then
             xf_iprn(i) = xf_iprn(i) + 500
         end if
         if( xf_svtype(i).eq.'I' ) then
             xf_iprn(i) = xf_iprn(i) + 600
         end if
      end do
      if ( debug_start.eq.-2 ) then 
         write(*,400) date, sectag, xf_msat,
     .               (xf_svtype(i),xf_iprn(i),i=1,xf_msat)
 400     format('RAW ',5i3, F6.1, 1x, i3, ' : ',100(1x,a1,1x,i3))
      endif     
 
****  Thats all
      return
      end


CTITLE READ_RINEX_OBS_2
 
      subroutine read_rinex_obs_2(eof, ep)

      implicit none
 
*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include '../includes/xfile_def.h'
      include 'track_com.h'

* PASSED
      integer*4 ep  ! Epoch number 

*   date(5)     - Ymdhm of observation
*   flags(10)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter 
*   jlk         - Dummy lock flag for reading SNR/Doppler entries
*   nl   -- Number of lines to be read 
      integer*4  flags(xf_maxsat), ierr, jerr, i,j,k, jlk
      integer*4  itmp, nl, trimlen

      integer*4  sys   ! GSNSS system number: 1 GPS, 2 Glonass, 3 Galileo, 4 Beidou
      
*   vals(10)     - Range values read
 
      real*8  vals(xf_maxsat)
 
*   eof     - Indicates end of file.
*   OK_data(xf_maxsat) -- Array set true if data on PRN seems
*     OK, otherwize set to .false. and is deleted from data
*     set.
 
      logical eof, OK_data(xf_maxsat)
      logical baddata(max_prn)   ! Flags reason for bad data
      logical proc_sys      ! Set true if we are processing this GNSS

      character*256 line
      character*1 cr
      
      cr = char(13)
 
 
*     Now loop over the data records.
      itmp = min(xf_ndat,5)
      nl = (xf_ndat-1)/5+1
      do i = 1, xf_msat
* MOD TAH 090114: Updated to allow for multiple lines
          do k = 1, nl 
             read(obs_lu,'(a)',iostat=ierr) line
             call report_error('IOSTAT',ierr,'read',obs_file,0,
     .                          'RX_RINEX_OBS_2')

             call sub_char(line,cr," ") 
             if( k.lt.nl ) then
                 itmp = min(5*k,xf_ndat)
             else
                 itmp = xf_ndat
             endif
             if( ierr.eq.0 )
     .       read(line,120, iostat=jerr) (vals(j), xf_ilck(j,i),
     .                flags(j), j =(k-1)*5+1,itmp)
 120         format(5(f14.3,i1,i1))
             call report_error('IOSTAT',jerr,'read',line,0,
     .                          'RX_RINEX_OBS_2')
          end do

          sys = int(xf_iprn(i)/100)+1

          do k =1, xf_mdat
             xf_obsv(k,i) = 0.0
             xf_isnr(k,i) = 0
* MOD TAH 180322: Updated from Multi-GNSS. Set for GNSS Type
             if( xf_dattyp(k,sys).gt.0 ) then
                 j = xf_dattyp(k,sys)
                 if( vals(j).ne.0.d0 ) then
                    xf_obsv(k,i) = vals(j)
                    xf_isnr(k,i) = flags(j)
                 end if
             end if
          end do 
                          
C            do j = 1 ,xf_ndat 
C               if( xf_dattyp(j).eq.k ) then
*                   Only save the value the input is non-zero
C                   if( vals(j).ne.0.d0 ) then
C                      xf_obsv(k,i) = vals(j)
C                      xf_isnr(k,i) = flags(j)
C                   end if
C               end if
C            end do
C         end do

* MOD TAH 990110: Now check that the data look basically
*         OK.  Checks are ranges and/or phases are zero.
          OK_data(i) = .true.  
* MOD TAH 100717: Check phase and then range
          baddata = .false.
          do k = 1, 2
             if( xf_obsv(k,i).eq.0.d0 ) OK_data(i) = .false.
             if( xf_obsv(k,i).eq.0.d0 ) baddata(1) = .true.
          end do
* MOD TAH 100717: Check P1 and C1
          if( xf_obsv(3,i).eq.0.d0 .and. xf_obsv(5,i).eq.0.d0 ) then
             OK_data(i) = .false.
             baddata(2) = .true.
          end if 
* MOD TAH 1404010: Add test in for C2 data (type 4: P2 type 6: C2)
          if( xf_obsv(4,i).eq.0.d0 .and. xf_obsv(6,i).eq.0.d0 ) then 
             OK_data(i) = .false.
             baddata(3) = .true.
          end if
!         if( baddata(3) ) 
!    .        print *,'CODES Num ',xf_ndat, ' TYPES ', xf_dattyp
* MOD TAH 100717: Check P1/C1 : These tests should be updated for C2
          if( xf_obsv(3,i).ne.0 .and. xf_obsv(4,i).ne.0 .and.
     .       abs(xf_obsv(3,i)-xf_obsv(4,i)).gt.200.d0  )
     .                                  OK_data(i) = .false.
          if( xf_obsv(3,i).ne.0 .and. xf_obsv(4,i).ne.0 .and.
     .       abs(xf_obsv(3,i)-xf_obsv(4,i)).gt.200.d0  )
     .                                  baddata(4) = .true.
          if( xf_obsv(5,i).ne.0 .and. xf_obsv(4,i).ne.0 .and.
     .       abs(xf_obsv(5,i)-xf_obsv(4,i)).gt.200.d0  )
     .                                  OK_data(i) = .false.
          if( xf_obsv(5,i).ne.0 .and. xf_obsv(4,i).ne.0 .and.
     .       abs(xf_obsv(5,i)-xf_obsv(4,i)).gt.200.d0  )
     .                                  baddata(5) = .true.

* MOD TAH 100719: Make sure L2 is not L1*60/77
          if( abs(xf_obsv(2,i)-xf_obsv(1,i)*60.d0/77.d0).lt.0.002 .and.
     .        xf_obsv(1,i).ne. 0.d0 ) then
              OK_data(i) = .false.
!              print *,'L2=L1*60/77 PRN ',xf_iprn(i),' EP ',ep,
!     .             ' Site ',site_names(obs_lu-50) 
              baddata(6) = .true.
          endif 

* MOD TAH 030729: See if in exclude list
          do k = 1, num_exclude
             if( xf_iprn(i).eq.ss_exclude(k) ) OK_data(i) = .false.
          end do 

* MOD TAH 040428: See if non-GPS (only output debg if GPS).
* MOM TAH 180322: See if GNSS system being prococessed 
          proc_sys = .true.
          call check_gnss( xf_iprn(i), tr_gnss, proc_sys )
          OK_data(i) = OK_data(i) .and. proc_sys 

          if( ep.ge.debug_start .and. ep.le. debug_end ) then
             write(*,210) ep, i, xf_iprn(i), proc_sys, OK_data(i), 
     .                tr_gnss
 210         format('CHECK_GNSS EP ',I5,' CH/PRN ',2I4,' STATE ',2L2,
     .              ' TR_GNSS ',a)
          end if

          if( .not.OK_data(i) .and. debug_start.eq.-2 
     .        .and. proc_sys ) then
             print *,'Err dP12 EP CH ',ep, i, 
     .           abs(xf_obsv(3,i)-xf_obsv(4,i)),
     .           xf_obsv(3,i),' dC1P2 ',
     .           abs(xf_obsv(5,i)-xf_obsv(4,i)),  xf_obsv(5,i),
     .           ' dL1L2 ',
     .           abs(xf_obsv(2,i)-xf_obsv(1,i)*60.d0/77.d0)
             print *,'BADDATA Logical ',ep, ' CH ',i,' BD ',
     .           baddata(1:8)
          end if
          

      end do

* MOD TAH 990110: Remove any data that has been flagged bad.
      i = 0
      do while (i.lt.xf_msat )
         i = i + 1
         if( .not.OK_data(i) ) then

*            This is bad measurement, so remove by moving
*            all data down and reducing xf_msat
             if( (ep.ge.debug_start .and. ep.le.debug_end .and.
     .           ep.ne.0).or. debug_start.eq.-2  ) then
                 write(*,220) i, xf_iprn(i), ep, (xf_obsv(k,i),  
     .                        xf_isnr(k,i), k=1, xf_ndat)
 220             format('Deleting CH ',I2,' PRN ',i3.2,' EP ',i6,
     .                 ' Data ',100(E13.3,1x,i4,1x))
             endif
*
             do j = i+1, xf_msat
                OK_data(j-1) = OK_data(j)
                xf_iprn(j-1) = xf_iprn(j)
                do k = 1, xf_mdat
                   xf_obsv(k,j-1) = xf_obsv(k,j) 
                   xf_isnr(k,j-1) = xf_isnr(k,j)
                end do
             end do
*            Now decrement the number of channels and the current
*            channel number (since we have just replaced this with
*            entry above).
             xf_msat = xf_msat - 1
             i = i - 1
         end if
      end do
      if ( debug_start.eq.-2 ) then
         write(*,290) ep, xf_msat
 290     format('DEBUG: EP ',i5,' MM_MSAT ',i3)
      endif

      if( ierr.ne.0 ) then
           eof = .true.
           call report_error('IOSTAT',ierr,'read',line,0,
     .                          'RX_RINEX_OBS_2')
      endif

      if( ep.ge.debug_start .and. ep.le.debug_end .and.
     .           ep.ne.0 ) then
*         DEBUG: Only first 8 values defined for GPS
          do i = 1, xf_msat
            write(*,300) ep, i, xf_iprn(i), xf_obsv(1:6,i)
 300        format('FIN EP ',i5,' CH ',i2,' PRN ',I3,' DATA ',
     .             10(F15.3,1x))
          end do
      end if 
      
****  Thats all

      return
      end
 
CTITLE CLEAN_CHAR
 
      subroutine clean_char(line, char_1, length)
 
      implicit none

*     This routine will clean the character in a line
      character*(*) line
      character*1 char_1
      integer*4 length, i,trimlen
      
      do i=1, trimlen(line)
        if (line(i:i).eq.char_1) line(i:i)=' '
      enddo
      
      return
      end
      
      
