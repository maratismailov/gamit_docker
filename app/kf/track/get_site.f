ctitle
 
      subroutine get_site(buffer,sit_num,itype,indx)

      implicit none
  
c
c     routine to extract the information form a "station" entry
c     in the Markov command file

c
c Variables
c ---------
c buffer  -- this is the line which has been read from the control
c        file with the station name extracted from the beginning of
c        of the line.
c sit_num -- the number of the site as given by the position of the
c        station name in the site_names array.
c itype   -- this integer tells the type of information to be read
c        from the line. (see subroutine assign_control.)
c
c Local variables
c ---------------
c label   -- integer variable which is used to convert itype to
c        executables statements through a computed goto.
c vals  -- real*8 temporary array which is used to store vals
c        from the line before their vals are converted to internal
c        units.
c
      include 'makexk_cmds_bd.h'
*                                          ! the parameter file
      include '../includes/xfile_def.h'
      include 'track_com.h'
c
      character*(*) buffer
 
      integer*4 sit_num, itype, label
 
c
*   err     - Error returned from multiread
*   indx   - counter used by multiread to keep track
*           - of position in command line.
*   jerr - Error checking if number
*   jndx - Pointer in string

      integer*4 err, indx,i, j, jerr, jndx
      integer*4 date(5), trimlen
      real*8 sec
 
*   vals(10)   - Temporary storage of vals read from command
*               - line.
 
      real*4 vals(10)
      real*8 vals8(10)
 
*   cvalue      - Dummy character string used in multiread
 
      character*256 cvalue
      character*2  code  ! DCB code (N/C/P) read from rcv_type line

      character site_nam*16, flg1*1, flg2*1
      integer*4 dlat,mlat, dlon, mlon
      real*8 slat,slon,radius

      integer*4 num   ! Number  of arguments passed in site_apr line
     
c.... compute the goto label based on the itype value.
      label = itype - site_start + 1

c.... goto the appropriate code
      goto(  100, 200, 300, 400, 500, 600, 700 ) label
c
c.... get the atmospheric delay offset
  100 continue
         call read_line(buffer,indx,'R8',err,vals8,cvalue)
         atm_offset(sit_num) = vals8(1)
         return

c.... atmosphere apriori vals
  200 continue

*        Get the apriori sigma
         call read_line(buffer,indx,'R4',err,vals,cvalue)
         apr_atm(sit_num) = vals(1)**2

*        Now get the process noise: This is m**2/epoch although
*        entry in m/(sqrt epoch)
         call read_line(buffer,indx,'R4',err,vals,cvalue)
         if( err.eq.0 ) mar_atm(sit_num) = vals(1)**2

*        Now get the height rate process noise.  This if for
*        aircraft that are radiply changing height
* MOD TAH 1.26: Added SCALE as option here instead of height rate
         jndx = indx
         call read_line(buffer,jndx,'CH',err,vals,cvalue) 
         call casefold(cvalue)
         if( cvalue(1:1).eq.'S' ) then  ! Version 1.26 Scale option
                                        ! selected.
             atm_scale(sit_num) = .true.
             mar_atm_hgt(sit_num) = 0.d0
         else        ! See if height rate passed
             call read_line(buffer,indx,'R4',err,vals,cvalue)
             if( err.eq.0 ) mar_atm_hgt(sit_num) = vals(1)**2
             atm_scale(sit_num) = .false.
         endif

      return

c.... site position aprioris
  300 continue
         num = -7
         call multiread(buffer,indx,'R8',err,vals8,cvalue,num)
c
c....    compute variances
	 do i = 1, 3
             site_int(i,sit_num)= vals8(i)
             site_apr(i,sit_num)= vals8(i)
         enddo
         if( num.eq.7 ) then
             call decyrs_to_mjd( vals8(7),site_ep(sit_num))
             do i = 1,3
                site_vel(i,sit_num) = vals8(i+3)
             end do
         else
             site_ep(sit_num) = 51544.0d0 ! Jan 1, 2000
             site_vel(:,sit_num) = 0.d0
         endif

         read_site_apr(sit_num)=.true.
       
      return

c TIMEDEP_PROCNS
c  Site    Sig XYZ (m/sqrt(t))  Start YY MM DD MN Sec End YY MM DD MN Sec
  400 continue
*        Increment counter
         num_tdep = num_tdep + 1
         j = num_tdep
         if( num_tdep.gt.5*max_site ) then
            call report_stat('FATAL','track','get_site',
     .         ' ','Too many timedep_procns entries',5*max_site)
         end if

*        Save the site number
         if( sit_num.gt.999 ) then
            nst_tdep(j) = -1
         else
            nst_tdep(j) = sit_num
         end if
****     Decode values
         call multiread(buffer,indx,'R4',err,vals,cvalue,3)
         do i = 1,3
            mar_tdep(i,j) = vals(i)**2
         end do
*        Start time
         call multiread(buffer,indx,'I4',err,date,cvalue,5)
         call multiread(buffer,indx,'R8',err,sec,cvalue,1)
         call ymdhms_to_mjd(date, sec, mjd_tdep(1,j))
*        End time
         call multiread(buffer,indx,'I4',err,date,cvalue,5)
         call multiread(buffer,indx,'R8',err,sec,cvalue,1)
         call ymdhms_to_mjd(date, sec, mjd_tdep(2,j))

C OLD NOT USED code
C         read(buffer, 410)   site_nam, flg1, dlat, mlat, slat,
C    .          flg2, dlon, mlon, slon, radius
C    
C         write(xf_posflg, '(a1,a1,"U")') flg1, flg2
C         call geoxyz(dlat,mlat,slat,flg1,dlon,mlon,slon,
C    .             flg2,radius,xf_pos(1), xf_pos(2), xf_pos(3))
C         write(xf_siteinfo, 420) site_nam, flg1, dlat, mlat, slat,
C    .         flg2, dlon, mlon, slon,radius, xf_rcvrsw,xf_swver    
C410      format(4x,a4,1x,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,
C    .         1x,f8.5,f13.4,2x,a3,1x,f5.2 )
C420      format(1x,a16,1x,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,
C    .         1x,f8.5,f13.4,2x,a3,1x,f5.2 )
C         do i = 1, 3
C             site_apr(i,sit_num)= xf_pos(i)
C         enddo
C         read_site_apr(sit_num)=.true.
      return
 
c.... site position noise level.  This sets the apriori sigma and 
*     process noise.  Units are m**2 and m**2/epoch.  NOTE: These
*     values are ignored for the F fixed sites.
* MOD TAH 100520: Added feature to constrain aposterori sigma
*     but putting POST after apriori sigma entries.
  500 continue
         call multiread(buffer,indx,'R4',err,vals,cvalue,3)
*        See if POST string is added
         jndx = indx
         call read_line(buffer,jndx,'CH',err, vals, cvalue)
         call casefold(cvalue) 
         if( cvalue(1:4).ne.'POST' ) then
            do i = 1,3
               apr_site(i,sit_num) = vals(i)**2
            end do

*           Now get the markov process noises
            call multiread(buffer,indx,'R4',err,vals,cvalue,3)
            do i = 1,3
               mar_site(i,sit_num) = vals(i)**2
            end do
         else
            do i = 1,3
               pst_site(i,sit_num) = vals(i)**2
            end do
         end if
      return 

c.... RCV_TYPE 
c          <Name> <code>
  600 continue

***      Allows specification of DCB type code for receivers at sites
         call GetWord(buffer,code, indx)
         call casefold(code)
         
         if( code(1:1).eq.'N' .or. code(1:1).eq.'C' .or.
     .       code(1:1).eq.'P' ) then
             rcv_type(sit_num) = code
         else
             write(*,620) code, buffer(1:trimlen(buffer))
 620         format('**WARNING** Unknown DCB code ',a1,
     .              ' RCV_TYPE line ',a)
         end if

      return
 

c.... get the antenna offset for ARP and L1 L2. 
* MOD TAH 061228: Change command structure to allow antenna name after
*     ARP values 
  700 continue
         call multiread(buffer,indx,'R8',err,vals8,
     .                   cvalue,3)
         do i = 1, 3
            site_offarp(i,sit_num)= vals8(i)
         end do

****     Get next value from line and see if number (old format)
         jndx = indx
         call read_line(buffer,jndx,'CH',err,vals8,cvalue)
         call check_num(cvalue,jerr)
         if( jerr .eq.0 ) then    ! Value is number, read rest of line
            call multiread(buffer,indx,'R8',err,vals8,
     .                      cvalue,6)
            do i = 1, 3
               sit_L12(i,1,sit_num)= vals8(i)
               sit_L12(i,2,sit_num)= vals8(i+3)
            enddo
         else   ! New format where antenna name string is passed (20-characters)
*           read entry under format control to allow for spaces in name
            ante_off = .true.
            cvalue = buffer(indx:)
            call casefold(cvalue)
            call trimlead(cvalue) 
            if( cvalue(17:20).eq.'    ' ) cvalue(17:20) = 'NONE'
            ant_name(sit_num) = cvalue
         end if  
****     See if receiver type passed as last argument
         indx = trimlen(buffer)-1
         call read_line(buffer,indx,'CH',err,vals8,rcv_type(sit_num))

      return

c.... This is the end of the routine
      end
 
c.........................................................................
