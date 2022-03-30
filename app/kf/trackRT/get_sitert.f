ctitle GET_SITERT
 
      subroutine get_sitert(buffer,sit_num,itype,indx)

      implicit none
 
      include '../includes/const_param.h' 
c
c     routine to extract the information form a "station" entry
c     in the Markov command file
* MOD TAH 091030: For use in TrackRT

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
      include 'trackRT_cmds.h'
*                                          ! the parameter file
      include 'trackRT.h'
c
      character*(*) buffer
 
      integer*4 sit_num, itype, label
 
c
*   err     - Error returned from multiread
*   indx   - counter used by multiread to keep track
*           - of position in command line.
*   jerr - Error checking if number
*   jndx - Pointer in string
*   num  - Number of values read from apr line

      integer*4 err, indx,i, jerr, jndx, num
      integer*4 trimlen
 
*   vals(10)   - Temporary storage of vals read from command
*               - line.
 
      real*4 vals(10)
      real*8 vals8(10)
      real*8 rot(3,3), loc(3)
 
*   cvalue      - Dummy character string used in multiread
 
      character*256 cvalue
      character*2  code  ! DCB code (N/C/P) read from rcv_type line

      character site_nam*16, flg1*1, flg2*1
      integer*4 dlat,mlat, dlon, mlon
      real*8 slat,slon,radius
     
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
         call read_line(buffer,indx,'R4',err,vals,cvalue)
         if( err.eq.0 ) mar_atm_hgt(sit_num) = vals(1)**2

        
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
*        compute init lat and long (used for CSV file header)
         call XYZ_to_geod(rot,site_int(1,sit_num), loc)
         init_geod(1,sit_num) = loc(1)*180/pi
         init_geod(2,sit_num) = loc(2)*180/pi
         init_geod(3,sit_num) = loc(3)
       
      return

c.... RCV_TYPE <Name> <code>
  400 continue

***      Allows specification of DCB type code for receivers at sites
         call GetWord(buffer,indx, code)
         call casefold(code)
         if( code(1:1).eq.'N' .or. code(1:1).eq.'C' .or.
     .       code(1:1).eq.'P' ) then
             rcv_type(sit_num) = code
         else
             write(*,420) code, buffer(1:trimlen(buffer))
 420         format('**WARNING** Unknown DCB code ',a1,
     .              ' RCV_TYPE line ',a)
         end if

      return
 
c.... site position noise level.  This sets the apriori sigma and 
*     process noise.  Units are m**2 and m**2/epoch.  NOTE: These
*     values are ignored for the F fixed sites.
  500 continue
         call multiread(buffer,indx,'R4',err,vals,cvalue,3)
         do i = 1,3
            apr_site(i,sit_num) = vals(i)**2
         end do

*        Now get the markov process noises
         call multiread(buffer,indx,'R4',err,vals,cvalue,3)
         do i = 1,3
            mar_site(i,sit_num) = vals(i)**2
         end do
      return 
c
c.... get the Pressure, temperature, Rh at the sites 
c     associated with their recorded height.  
  600 continue
*        Command no longer used.    
          write(*,'(a)') '**WARNING** MET_DATA command not implemented'
      return

c.... get the antenna offset for ARP and L1 L2. 
* ANTE_OFF 
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
 
