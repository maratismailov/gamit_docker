CTITLE GLOBK_SUM
 
      subroutine globk_sum ( lu_out )
 
      implicit none

 
*     Routine to summarize the data in the global solution.  This
*     used to give an idea of what is in the solution for input data.
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i,j,k       - Loop counters
*   ierr        - IOSTAT write error
*   lu_out      - Output Lu number
*   trimlen     - HP function for length of a string
 
      integer*4 i,j, ierr, lu_out, trimlen
 
*   dec             - Dec of source in main memory
*   ded, dem, des   - Degrees minutes seconds of declination
*   ra              - RA in millitimeseconds
*   rah, ram, ras   - Hours minutes seconds part of RA
 
 
      real*8 dec, ded, dem, des, ra, rah, ram, ras
 
***** Tell user about the solution.
 
      write(lu_out, 100, iostat=ierr ) list_file(1:trimlen(list_file))
  100 format(//' GLOBK: ',t40,'List ',a,/,
     .         ' ------', t40,'----'       )
 
      write(lu_out, 150, iostat=ierr) gnum_sites, gnum_sources,
     .                                num_glb_sol
  150 format(/' There are ',i3,' sites and ',i3,' sources in ',
     .        i5,' experiments')
 
      write(lu_out, 170, iostat=ierr) gnum_obs,
     .      (gdelete_count(i),i=1,max_edit_types)
  170 format(' There ',i6,' observations in the solution',/,
     .       ' with deletes ',10i6,/,
     .       '              ', 6i6)
 
***** List site names and postions
 
      write(lu_out, 200, iostat=ierr)
  200 format(/' SITES',/,'  Name ',t18,'X (m)',t35,'Y (m)',t52,
     .                   'Z (m)',t64,' Axis offset (m)')
 
      do i = 1, gnum_sites
          write(lu_out,250, iostat=ierr) gsite_names(i),
     .                     (apr_val_site(j,1,i),j=1,3),
     .                      apr_val_axo(1,i)
 250      format(1x,a8,t10,4(f16.4,1x))
      end do
 
***** List source names and positions
 
      write(lu_out, 300, iostat=ierr)
 300  format(/'SOURCES',/,'  Name',t17,'RA (hms)',t31,'Dec (dms)')
 
      do i = 1, gnum_sources
*                                             ! Convert to millitimesecs
          ra = apr_val_source(1,1,i)/15.d0
 
*                                             ! Convert to hrs, min,sec
          call mas_to_dms( ra,  rah, ram, ras)
 
          dec = apr_val_source(2,1,i)
          call mas_to_dms( dec, ded, dem, des )
 
          write(lu_out, 350, iostat= ierr) gsource_names(i),
     .            rah, ram, ras,  ded, dem, des
  350     format(1x,a8,t12,i3,1x,i2,1x,f9.6,1x,i3,1x,i2,1x,f9.6)
      end do
 
***** Thats all
      return
      end
 
