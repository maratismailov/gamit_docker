ctitle
 
      subroutine get_obsfile(buffer,itype)

      implicit none

c     Gang Chen, MIT, 96-97  
c
c     routine to extract the obs. data file name and information from 
c     the  parameter file
c     MOD, 09/97, G. Chen, add mark to indicate the type of input data
c              obs_file_type "R", "X", "C"        
      include 'makexk_cmds_bd.h'
      include 'track_com.h'
      
c Variables
c ---------
c buffer  -- this is the line which has been read from the control
c        file with the station name extracted from the beginning of
c        of the line.

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
c cvals -- character*20 temporary array which is used to store words

      character*(*) buffer
 
      integer*4 itype, label, trimlen
      
      character cvals(10)*40
      real*4 vals
c     
 
c
*   err     - Error returned from multiread
*   indx   - counter used by multiread to keep track
*           - of position in command line.
      integer*4 ierr, indx
      
c.... compute the goto label based on the itype value.
      label = itype - file_start + 1
*                    ! Start a beginning of line.
*
*     Read from the first column of line
      indx = 1
      if (trimlen(buffer).eq.0) then
        itype = -3
        return
      endif
      
c.... goto the appropriate code
      goto(  100, 200) label
            
c.... Decode the input-file line      
 100  continue

*        Get the observation file descriptors
         call multiread(buffer,indx,'CH',ierr,vals,cvals,3)
         if (ierr.eq.0) then
             num_site = num_site + 1
             if( num_site.gt.max_site ) then
                 call report_stat('FATAL','TRACK','get_obsfile',' ',
     .                            'Too many sites: Max allowed',
     .                            max_site)
             end if
             site_names(num_site) = cvals(1)
             obs_file(num_site) = cvals(2)
             obs_file_type(num_site) = "R"
             swvers(num_site) = '0.0'    
             call casefold(cvals(3))
             if (cvals(3).eq.'F') site_type(num_site) = 0
             if (cvals(3).eq.'K') site_type(num_site) = 1
         end if

      return

****  Option to input xfile type:
 200  continue
*       Option no longer supported
  
        return

****  Thats all
      end       

