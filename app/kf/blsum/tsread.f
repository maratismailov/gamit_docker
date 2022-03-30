      subroutine merge_ts(ns)

      implicit none

*     Routine to merge the new time series values with the old ones
*     If a time matches with an old value, the older value will be replaced.

      include 'tssum.h'

      integer*4 ns  ! code site number

      integer*4 i,j,k, n, l

      logical replaced ! Logical set true if entry is being replaced with
                       ! a new one.


      do i = 1, num_ent
         if( in_cs(i).eq.ns) then

****        The site code matches, scan the old time series values
*           to see we need to replace the value, insert a new one or
*           append to the old list (later would be most common).
*           Times may differ by a minute or so; so test with toleranace
C           if( num_ts.gt.0 .and. in_mjd(i).le.ts_mjd(num_ts) ) then
            replaced = .false.
            if( num_ts.gt.0 .and. 
     .          in_mjd(i)-ts_mjd(num_ts).lt.1.d-3 ) then
*               This entry is earlier.  Find out where to put it.
                n = num_ts
                do while ( n.ge.1 ) 
C                    if ( ts_mjd(n).eq.in_mjd(i) ) then
                    if ( abs(ts_mjd(n)-in_mjd(i)).le.1.d-3 ) then
*                      Matches the time, so replace the old values
*                      with the new
                       k = n
                       n = 0
                       replaced = .true.  
                    elseif ( ts_mjd(n).lt.in_mjd(i) .or.n.eq.1) then
*                      Need to move up all the old values and insert
*                      the new entry in at n+1
                       k = n+1
* MOD TAH 131203:      If inserted at the first point, set k = 1
                       if( n.eq.1 .and. ts_mjd(n).ge.in_mjd(i)) k = 1

*                      Move up the old values
                       do l = num_ts+1, k, -1
                          ts_mjd(l) = ts_mjd(l-1)
                          ts_type(l) = ts_type(l-1)
                          do j = 1, 3
                             ts_xyz(j,l) = ts_xyz(j,l-1)
                             ts_neu(j,l) = ts_neu(j,l-1)
                             ts_llu(j,l) = ts_llu(j,l-1)
                          enddo
                          do j = 1, 6
                             ts_xyz_std(j,l) = ts_xyz_std(j,l-1)
                             ts_neu_std(j,l) = ts_neu_std(j,l-1)
                          enddo
                       end do
                       num_ts = num_ts + 1
                       n = 0
                    else
                       n = n -1
                    end if

                end do
            else 
*               Normal case of just adding new results to end
                num_ts = num_ts + 1
                k = num_ts
            endif

*           Now add the new values
            ts_mjd(k) = in_mjd(i)

*           See which program and use this to set code
            if( tsprog(1:5).eq.'tscon' ) then
                ts_type(k) = in_type(i)
            else
*               Get which code we should add for this point.  
                if( replaced .and. ts_ref_type.eq.'suppl' ) then
                    ts_type(k) = 'suppf'
                elseif( replaced .and. ts_ref_type.eq.'supp6' ) then
                    ts_type(k) = 'sup6f'
                elseif( replaced .and. ts_ref_type.eq.'campd' ) then
                    ts_type(k) = 'campf'
                else
                    ts_type(k) = ts_ref_type
                endif
            endif

            do j = 1, 3
               ts_xyz(j,k) = in_xyz(j,i)
               ts_neu(j,k) = in_neu(j,i)
               ts_llu(j,k) = in_llu(j,i)
            enddo
            do j = 1, 6
               ts_xyz_std(j,k) = in_xyz_std(j,i)
               ts_neu_std(j,k) = in_neu_std(j,i)
            enddo
  
         end if
      end do
****  Thats all
      return
      end

CTITLE GEN_REFDATA

      subroutine gen_refdata(ns)

      implicit none

      include 'tssum.h'

      integer*4 ns  ! code site number

      real*8 av_xyz(3)
      integer*4 na, i, j

****  Routine to generate reference XYZ
      do j = 1,3 
        av_xyz(j) = 0.0d0
      end do
      na = 0

      do i = 1, num_ent
         if( in_cs(i).eq.ns) then
             if( in_xyz_std(1,i).lt.0.1d0 .and.
     .           in_xyz_std(2,i).lt.0.1d0 .and.
     .           in_xyz_std(3,i).lt.0.1d0 ) then 
                na = na + 1
                do j = 1,3
                   av_xyz(j) = av_xyz(j) + in_xyz(j,i)
                end do
             end if
         end if
      end do
      if( na.gt.0 ) then
          do j = 1,3
             ref_xyz(j) = av_xyz(j)/na
          end do
      else
*         Nothing with sigmas less than 10 cm: So what to do?
*         Use all available data?
          do i = 1, num_ent
             if( in_cs(i).eq.ns) then
                 na = na + 1
                 do j = 1,3
                    av_xyz(j) = av_xyz(j) + in_xyz(j,i)
                 end do
              end if
          end do
          if( na.gt.0 ) then
              do j = 1,3
                 ref_xyz(j) = av_xyz(j)/na
              end do
          end if
      endif
 

****  Thats all
      return
      end

CTITLE REMOVE_EJMP

      subroutine remove_ejmp(ns)

      implicit none

*     Subroutine to remove east jumps that can occur when a site
*     moves too far in latitude.  Usually only a problem for ice
*     sites.  Here we compare the reference lat with the latitude
*     in the individual entries

      include '../includes/const_param.h'
      include 'tssum.h'

      integer*4 ns  ! code site number; if site 0 is used, the ts_ values are
                    ! used.

      integer*4 i

      real*8 locref(3), locdat(3)  ! Geod Co-lat, long and heigt
                                   ! with reference latitude and 
                                   ! day value of latitude
      real*8 rot(3,3)
      real*8 neuref(3), neudat(3)  ! NEU values for reference and
                                   ! daily latiude.  If no latitude
                                   ! induced East jump, then should
                                   ! be the same
      real*8 ref_neu_base(3)       ! Reference position NEU values.
      real*8 dEast, old_dEast      ! East change and previous east change
                                   ! (latter used for output)

      logical debug
      debug = .false.
*     if( in_code(abs(ns)).eq.'P630' ) debug = .true.
      if( ns.lt.0 ) ns = 0         ! Reset back to ts values 
*
***** Get the reference positions for this site.
      call XYZ_to_GEOD(rot,ref_xyz, ref_llu) 
      call loc_to_geod(ref_llu, ref_neu_base)

*     Loop over entries finding the one for this site
      old_dEast = 0.0d0

      if( ns.gt. 0 ) then 
         do i = 1, num_ent
            if( in_cs(i).eq.ns) then
*              Test to see if change.  Call loc_to_geod with
*              actual and refence latitude
               locdat(1) = pi/2-in_llu(1,i)*pi/180
               locref(2) = in_llu(2,i)*pi/180
               locdat(2) = in_llu(2,i)*pi/180
               locref(3) = in_llu(3,i)
               locdat(3) = in_llu(3,i)
*              Now change reference lat in locref
               locref(1) = ref_llu(1)
*              Convert to NEU
               call loc_to_geod(locref, neuref)
               call loc_to_geod(locdat, neudat)
               if( debug ) then 
                   print *,' Ref ENT ',i, ' NEU ', neuref, in_mjd(i) 
                   print *,' Dat ENT ',i, ' NEU ', neudat, num_ent
               endif

*              Get the East difference
               dEast = neuref(2)-neudat(2)
               if( abs(dEast).gt.1.d-4 ) then
! MOD TAH 131101: Reset the NEU value to be consistent.
!                 in_neu(2,i) = in_neu(2,i)+dEast
!                 write(*,215) i, neuref(2), neudat(2), 
!     .                ref_neu_base(2), in_neu(2,i)
! 215             format('EP ',i4,' Easts ',3F15.4, ' IN ',G15.4)
                   in_neu(2,i) = neuref(2)
!                  ts_neu(2,i) = in_neu(2,i)
               endif
               if( abs(dEast-old_dEast).gt.1.d-2) then
                  if( debug )
     .            write(*,220) in_code(ns), in_mjd(i), dEast
 220              format('East Offset at ',a4,' MJD ',F12.4,
     .                   ' dEast ',F15.3,' m')
                  ts_code = in_code(ns)   ! Save name 
               endif
               old_dEast = dEast
            endif
          end do
       else      ! Call from tsfit where the ts_values are being
                 ! used

         do i = 1, num_ts
*            Test to see if change.  Call loc_to_geod with
*            actual and refence latitude
             locdat(1) = pi/2-ts_llu(1,i)*pi/180
             locref(2) = ts_llu(2,i)*pi/180
             locdat(2) = ts_llu(2,i)*pi/180
             locref(3) = ts_llu(3,i)
             locdat(3) = ts_llu(3,i)
*            Now change reference lat in locref
             locref(1) = ref_llu(1)
*            Convert to NEU
             call loc_to_geod(locref, neuref)
             call loc_to_geod(locdat, neudat)
             if( debug ) then 
                  print *,' Ref TSF ',i, ' NEU ', neuref, ts_mjd(i) 
                  print *,' Dat TSF ',i, ' NEU ', neudat, num_ts
              endif

*            Get the East difference
             dEast = neuref(2)-neudat(2)
             if( abs(dEast).gt.1.d-4 ) then
!               ts_neu(2,i) = ts_neu(2,i)+dEast
! MOD TAH 131101: Reset the NEU value to be consistent.
                ts_neu(2,i) = neuref(2)
             endif
             if( abs(dEast-old_dEast).gt.1.d-2) then
                if( debug ) 
     .          write(*,230) ts_code, ts_mjd(i), dEast
 230            format('East Offset at ',a4,' MJD ',F12.4,
     .                 ' dEast ',F15.3,' m: TSFIT')
             endif
             old_dEast = dEast
            
          end do
       end if

*****  Thats all
       return 
       end

CTITLE SAVE_TS

      subroutine save_ts(in,out)

      implicit none

*     Routine to save time series values from index in to index out
*     If out is zero, save to or from save arrrays

      include 'tssum.h'

      integer*4 in, out

* LOCAL
      integer*4 i


****  See which way we are saving
      if( in.gt.0 .and. out.gt.0 ) then
          ts_mjd(out) = ts_mjd(in)
          ts_type(out) = ts_type(in)
          do i = 1,3
             ts_xyz(i,out) = ts_xyz(i,in)
             ts_llu(i,out) = ts_llu(i,in)
             ts_neu(i,out) = ts_neu(i,in)
          enddo
          do i = 1,6
             ts_xyz_std(i,out) = ts_xyz_std(i,in)
             ts_neu_std(i,out) = ts_neu_std(i,in)
          enddo
      elseif( out .eq. 0 ) then
          sv_mjd = ts_mjd(in)
          sv_type = ts_type(in)
          do i = 1,3
             sv_xyz(i) = ts_xyz(i,in)
             sv_llu(i) = ts_llu(i,in)
             sv_neu(i) = ts_neu(i,in)
          enddo
          do i = 1,6
             sv_xyz_std(i) = ts_xyz_std(i,in)
             sv_neu_std(i) = ts_neu_std(i,in)
          enddo
      elseif( in.eq.0 ) then
          ts_mjd(out) = sv_mjd
          ts_type(out) = sv_type
          do i = 1,3
             ts_xyz(i,out) = sv_xyz(i)
             ts_llu(i,out) = sv_llu(i)
             ts_neu(i,out) = sv_neu(i)
          enddo
          do i = 1,6
             ts_xyz_std(i,out) = sv_xyz_std(i)
             ts_neu_std(i,out) = sv_neu_std(i)
          enddo
      else
          call report_error('PROG',-1,'sav','timeseries',0,'save_ts')
      endif
      return
      end

