CTITLE REPORT_NONSEC

      subroutine report_nonsec( iout )

      implicit none 

*     Routine to report the non-secular terms used in the globk
*     analysis

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED VARIABLES
* iout  -- Output unit number

      integer*4 iout

* LOCAL VARIABLES
* it -- Loop over types of non-secular terms
* is -- Site number
* i, j -- Loop counters
* date(5) -- YMDMN of non-secular term

      integer*4 it, is, i,j, date(5)

* dsol_xyz(3), dsol_neu(3) -- Changes in position at epoch for the
*    non-secular term
* loc_coord(3) -- Lat,long, and height
* rot_mat(3,3) -- Rotation matrix
* vals_neu(6)  -- Values of coefficients in NE and U
* sectag       -- Seconds tag for date conversion

      real*8 dsol_xyz(3), dsol_neu(3), loc_coord(3), rot_mat(3,3),
     .       vals_neu(6), sectag
      real*8 deps  !  Small change in time to ensure minute rounds up.

* header -- Set true once header is written
      logical header

      deps = 1.d-5   !  ~ 1sec.

****  OK, Loop over all the types and report them in blocks
      if( num_nonsec.eq.0 ) RETURN

      do it = 1, 4
         header = .false.
         do i = 1, num_nonsec

*           See if type is correct
            if( param_nonsec(2,i).eq.it .and. it.eq. 1 ) then
*               OK, we have an entry.  If we have not write
*               header write it now.
                if( .not.header ) then
                   write(iout,120) nonsec_types(it)
 120               format(/,' NON-SECULAR ',a8,' REPORT',/,
     .                    '          Site     Type      ',
     .                    'YY  MM DD HR MN     Not      ',
     .                    '   North              East   ',
     .                    '              Height         ',
     .                    'N      E      U',/,48x,
     .                    'Used    Offset     Rate    ',
     .                    'Offset     Rate    Offset  ',
     .                    '   Rate          at epoch (m)' ) 
                   header = .true.
                end if

*               Compute the contribution at the epoch of the soluiton
                is = param_nonsec(1,i)
                call eval_nonsec(is, gepoch_out, 1, param_nonsec(1,i),
     .                 apr_val_nonsec(1,i), dsol_xyz,0)
                call rotate_geod(dsol_xyz, dsol_neu,'XYZ','NEU',
     .                 apr_val_site(1,1,is), loc_coord, rot_mat)
*               Convert the XYZ values back to NEU for output.
                call nonsec_convert('TONEU',2,apr_val_nonsec(3,i),
     .                 vals_neu,apr_val_site(1,1,is))
                call jd_to_ymdhms(apr_val_nonsec(1,i)+deps,date,sectag)
                if( dsol_neu(1)**2+dsol_neu(2)**2+dsol_neu(3)**2 .gt.
     .              1.d-12 )
     .          write(iout,160) gsite_names(is), nonsec_types(it),
     .               date, apr_val_nonsec(2,i), vals_neu, dsol_neu
 160            format(' EXTENDED ',a8,1x,a8,1x,i5,4i3,F8.3,6F9.5,1x,
     .                 3F9.5)

             else if( param_nonsec(2,i).eq.it .and. it.eq.2 ) then
*               OK, we have an entry.  If we have not write
*               header write it now. 
*               PERODIC TERMS
                if( .not.header ) then
                   write(iout,220) nonsec_types(it)
 220               format(/,' NON-SECULAR ',a8,' REPORT',/,
     .               11x,'Site     Type     YY  MM DD HR MN',
     .                   '  Period         North              East',
     .                   '                 Height         ',
     .                   'N      E      U',/, 46x,
     .                   '(days)      Cos     Sin        Cos     Sin',
     .                   '          Cos    Sin          at epoch (m)')
                   header = .true.
                end if

*               Compute the contribution at the epoch of the soluiton
                is = param_nonsec(1,i)
                call eval_nonsec(is, gepoch_out, 1, param_nonsec(1,i),
     .                 apr_val_nonsec(1,i), dsol_xyz,0)
                call rotate_geod(dsol_xyz, dsol_neu,'XYZ','NEU',
     .                 apr_val_site(1,1,is), loc_coord, rot_mat)
*               Convert the XYZ values back to NEU for output.
                call nonsec_convert('TONEU',2,apr_val_nonsec(3,i),
     .                 vals_neu,apr_val_site(1,1,is))
                call jd_to_ymdhms(apr_val_nonsec(1,i)+deps,date,sectag)
                if( dsol_neu(1)**2+dsol_neu(2)**2+dsol_neu(3)**2 .gt.
     .              1.d-12 )
     .          write(iout,260) gsite_names(is), nonsec_types(it),
     .               date, apr_val_nonsec(2,i), vals_neu, dsol_neu
 260            format(' EXTENDED ',a8,1x,a8,1x,i5,4i3,F8.3,6F9.5,1x,
     .                 3F9.5)

             else if( param_nonsec(2,i).eq.it .and. it.eq.3 ) then
*               OK, we have an entry.  If we have not write
*               header write it now. 
*               EXPONENTIAL TERMS
                if( .not.header ) then
                   write(iout,320) nonsec_types(it)
 320               format(/,' NON-SECULAR ',a8,' REPORT',/,
     .                11x,'Site     Type     YY  MM DD HR MN',
     .                    '       Decay T   North     East   Height',
     .                    '       N       E        U',/, 46x,
     .                    '        (days)                           ',
     .                    '           at epoch (m)')

                   header = .true.
                end if

*               Compute the contribution at the epoch of the soluiton
                is = param_nonsec(1,i)
                call eval_nonsec(is, gepoch_out, 1, param_nonsec(1,i),
     .                 apr_val_nonsec(1,i), dsol_xyz,0)
                call rotate_geod(dsol_xyz, dsol_neu,'XYZ','NEU',
     .                 apr_val_site(1,1,is), loc_coord, rot_mat)
*               Convert the XYZ values back to NEU for output.
                call nonsec_convert('TONEU',1,apr_val_nonsec(3,i),
     .                 vals_neu,apr_val_site(1,1,is))
                call jd_to_ymdhms(apr_val_nonsec(1,i)+deps,date,sectag)
                if( dsol_neu(1)**2+dsol_neu(2)**2+dsol_neu(3)**2 .gt.
     .              1.d-12 )
     .          write(iout,360) gsite_names(is), nonsec_types(it),
     .               date, apr_val_nonsec(2,i), (vals_neu(j),j=1,3),
     .                dsol_neu
 360            format(' EXTENDED ',a8,1x,a8,1x,i5,4i3,F12.3,3F9.5,1x,
     .                 3F9.5)
             else if( param_nonsec(2,i).eq.it .and. it.eq.4 ) then
*               OK, we have an entry.  If we have not write
*               header write it now. 
*               LOGARITHMIC TERMS
                if( .not.header ) then
                   write(iout,420) nonsec_types(it)
 420               format(/,' NON-SECULAR ',a8,' REPORT',/,
     .                11x,'Site     Type     YY  MM DD HR MN',
     .                    '       Decay T   North     East   Height',
     .                    '       N       E        U',/, 46x,
     .                    '        (days)                           ',
     .                    '           at epoch (m)')

                   header = .true.
                end if

*               Compute the contribution at the epoch of the soluiton
                is = param_nonsec(1,i)
                call eval_nonsec(is, gepoch_out, 1, param_nonsec(1,i),
     .                 apr_val_nonsec(1,i), dsol_xyz,0)
                call rotate_geod(dsol_xyz, dsol_neu,'XYZ','NEU',
     .                 apr_val_site(1,1,is), loc_coord, rot_mat)
*               Convert the XYZ values back to NEU for output.
                call nonsec_convert('TONEU',1,apr_val_nonsec(3,i),
     .                 vals_neu,apr_val_site(1,1,is))
                call jd_to_ymdhms(apr_val_nonsec(1,i)+deps,date,sectag)
                if( dsol_neu(1)**2+dsol_neu(2)**2+dsol_neu(3)**2 .gt.
     .              1.d-12 )
     .          write(iout,460) gsite_names(is), nonsec_types(it),
     .               date, apr_val_nonsec(2,i), (vals_neu(j),j=1,3),
     .                dsol_neu
 460            format(' EXTENDED ',a8,1x,a8,1x,i5,4i3,F12.3,3F9.5,1x,
     .                 3F9.5)

             end if
          end do
      end do

****  Thats all
      return
      end



