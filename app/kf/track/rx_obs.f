CTITLE RX_OBS
 
      subroutine read_xf_obs_1(eof,buffer, in_code)

      implicit none
 
*     This routine will read the next group of ranges from the x file
*     eof if returned true if we run out of data.
*     Split this to 2 parts, 
*            read_rinex_obs_1 : read in date information
*            read_rinex_obs_2 : read in obs. data
*     Gang Chen  971011
 
      include '../includes/const_param.h'
      include '../includes/xfile_def.h'

*   in_code     -  0 normal read from file; 
*		-  1  from self buffer
      integer*4 in_code
 
*   date(5)     - hour minute of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error

      integer*4  date(5), ierr,  trimlen
 

      integer*4 hr, min, dlat,mlat, dlon, mlon
      real*8 slat,slon,height,sec, rot_matrix(3,3)

       
*   eof     - Indicates end of file.
 
      logical eof
 
      character*256 line, buffer
*      character*1 xf_outflg(1:1),xf_outflg(2:2)

      eof =.false.

      if(in_code.eq.0) then 
****  read the empty line
         read(obs_lu,'(a)', iostat=ierr) line

****  Read in the next line
         read(obs_lu,'(a)', iostat=ierr) line
         write(buffer,'(a256)') line     
         if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
             write(*,*) " here ",ierr,trimlen(line)
             eof = .true.
             RETURN
         end if
      else

****  Read in the next line from buffer
          read(buffer,'(a256)') line 
          in_code = 0    
c          write(*,*) "line",line     
      endif
      
 
*     check if empty record
      read(line,*,iostat=ierr) xf_iepoch, xf_msat
c
      if (xf_msat.eq. 0) then
        return
      end if
      
*     Decode the line
      read(line,100,iostat=ierr) xf_iepoch, xf_msat, xf_doy , 
     .    hr, min, sec,xf_kflag, xf_sitnam,
     .    xf_outflg(1:1), dlat, mlat, slat,
     .    xf_outflg(2:2), dlon, mlon, slon, height,
     .    xf_offsL1,xf_offsL2

 100  format(2I4,I5,I4,2I3,F11.7,1x,i2,1x,a4,1x,
     .          a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4,
     .          1x,3F8.4,3X,3F8.4)
       xf_sod = hr*3600 + min*60 + sec
       date(1) = xf_doy(1)
       date(2) = 1
       date(3) = xf_doy(2)
       date(4) = hr
       date(5) = min
       call ymdhms_to_jd( date, sec, xf_jd)
       call ymdhms_to_mjd( date, sec, xf_mjd )
       
c have problem with yds_to_jd
c       call ymdhms_to_jd( date, sec, xjd2)
c       call yds_to_jd(xf_doy(1), xf_doy(2), xf_sod, xf_jd) 
c       write(*,*) "yr", xf_doy, xf_sod, xf_jd   , xjd1,xjd2  

       xf_outflg(3:3) = 'U'

c       call geoxyz(dlat,mlat,slat,xf_outflg(1:1),dlon,mlon,slon,
c     .          xf_outflg(2:2),height, xf_pos(1), xf_pos(2), xf_pos(3))
       xf_posflg='XYZ'
          call XYZ_to_NEU(rot_matrix,xf_pos,xf_kneu)
c          write(*,*) xf_pos, xf_kneu, height

      if( ierr.ne.0 ) eof = .true.
 
****  Thats all
      return
      end


CTITLE READ_X_OBS
 
      subroutine read_xf_obs_2(eof)
 
      implicit none

*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include '../includes/xfile_def.h'

*   eof     - Indicates end of file.
 
      logical eof
 
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
      integer*4  ierr, i, k, ischan
 
      eof = .false.
      
      if (xf_msat.eq.0) return
      
*     Now loop over the data records.
      do i = 1, xf_msat
          read(obs_lu,120, iostat=ierr) xf_ierfl(i),ischan, 
     .        xf_obsv(1,i), xf_isnr(1,i),
     .        xf_obsv(2,i), xf_isnr(2,i),
     .        xf_obsv(3,i), xf_obsv(4,i)
          do k = 1, 2
             xf_obsv(k,i) =  - xf_obsv(k,i)
          end do  
          xf_iprn(i) = xf_prn(ischan)
c          write(*,*) i,xf_iprn(i), (xf_obsv(j,i), j=1,4)
             
      enddo
 120           format(10x,2i2,1x,d22.15,1x,i3,1x,d22.15,1x,i3,                 
     .            2x,d22.15,2x,d22.15) 


      if( ierr.ne.0 ) eof = .true.
 
****  Thats all
      return
      end
 
