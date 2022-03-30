      Program extend_iers

*     Program to read an IERS formatted polar motion/UT1 table and
*     to write globk and gamit polar motion/UT1 tables after extending
*     the 2-days forwards and backwards from the from the last and first
*     points.  This extension allows the tables to be used in gamit programs.
*
*     Runstring for program is:
*     % extend_iers <input file> <extent>
*     where <input file> is an IERS formatted PMU file (can be obtained
*           by grep 'IERS' <globk.prt/org file> | awk '{print substr($0,5)}'
*           <extent> is an extent to be put on the output files.  
*  
*     The output files are:
*     pmu.<root>  -- Globk PMU file
*     pole.<root> -- GAMIT pole table
*     ut1.<root>  -- GAMIT UT1 table. 

* PROGRAM PARAMETERS

* max_ent  -- Max number of lines in IERS tables

      integer*4 max_ent
      parameter ( max_ent = 10000 )

*
* PROGRAM VARIABLES      
* ----------------- 

* len_run, rcpar -- Length of runstring and function to return
* trimlen    -- Length of string
* ierr, jerr -- IOSTAT errors
* num_ent, nm -- number of entries found in IERS tables and short version
* ixpr(4), ispr(4)  -- Pole xy offset (1,2) and rates (3,4) and their
*               sigmas (as integers)
* iutr(2), isut(2)  -- UT1 and LOD and their sigmas
* num_sites, nf, num_svs -- Number of stations, fiducials and satellites
* i,j        -- Loop counters
* npr        -- Number of values per gamit table line
* tai_utc    -- Number of leap seconds
* date(5)    -- ymd h min for date format
* igxp(6), igyp(6) -- X and Y pole positions values to be output
* igut(6)    -- UT1 values to be output

      integer*4 len_run, rcpar, trimlen, ierr, jerr, num_ent, nm,
     .          ixpr(4), ispr(4), iutr(2), isut(2), num_sites, 
     .          nf, num_svs, i,j, npr, tai_utc,  date(5), 
     .          igxp(6), igyp(6), igut(6)

* rho(3)     -- Correlations between pole and ut1 values
* inf(3)     -- Information for interpolator: Start MJD, spacing, and
*               number of values
* pts(3,max_ent) -- Points to be interpolated (X, Y pole and UT1)
* rts(3,max_ent) -- Rates for points (only used is table < 4 entries 
*               long)
* jd, mjd   -- Generic JD and MJD
* djd       -- Difference in MJD from expected value
* bjd       -- MJD at start of each line
* st_mjd, en_mjd -- Start and stop MJD for gamit tables
* ut1_mult, pol_mult -- Multipliers for ut1 and pole values
* sec_tag   -- Seconds tag for epoch
* value(3,2) -- Interpolated values and their rates

      real*8 rho(3), inf(3), pts(3,max_ent), rts(3,max_ent),
     .       jd, mjd, bjd, st_mjd, en_mjd, ut1_mult, pol_mult, 
     .       sec_tag, value(3,2), djd
      integer*4 ibjd


* in_iers    -- Input IERS file name
* extent     -- Extent for naming output files
* out_globk  -- Output globk file name
* out_ut1_gamit -- Output name for gamit ut1 file
* out_pole_gamit -- Output name for gamit pole file
* line       -- Line read from file

      character*128 in_iers, extent, out_globk, out_ut1_gamit,
     .              out_pole_gamit, line

* ut1_format, pol_format -- Formats for UT1 and pole entries

      character*32 ut1_format
      character*35 pol_format 
* 

****  Get the runsting for the program
      len_run = rcpar(1,in_iers)
      if( len_run.eq.0 ) then
          call proper_runstring('extend_iers.hlp','extend_iers',1)
      else
          open(100,file=in_iers, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open',in_iers,1,
     .                      'extend_iers')
      end if

*     Get the extend string
      len_run = rcpar(2,extent)
      if( len_run.eq.0 ) then
          call proper_runstring('extend_iers.hlp','extend_iers',1)
      endif

      out_globk = 'pmu.' // extent
      out_pole_gamit = 'pole.' // extent
      out_ut1_gamit  = 'ut1.'  // extent
      write(*,100) out_globk(1:trimlen(out_globk)),
     .             out_ut1_gamit(1:trimlen(out_ut1_gamit)),
     .             out_pole_gamit(1:trimlen(out_pole_gamit)),
     .             in_iers(1:trimlen(in_iers))

100   format('Creating globk pmu file ',a,' and GAMIT files ',a,
     .       ' and ',a,' from ',a)

      open(200,file=out_globk, status='unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'creat',out_globk,1,
     .                  'extend_iers')
      open(201,file=out_ut1_gamit, status='unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'creat',out_ut1_gamit,1,
     .                  'extend_iers')
      open(202,file=out_pole_gamit, status='unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'creat',out_pole_gamit,1,
     .                  'extend_iers')

****  Read the first four entries in the file.  If we don't get four 
*     values then treat as special case and we will use the rates to
*     extrapolate the table.
      
      num_ent = 0
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr) line
*         Note: in the IERS format, the comment lines start with a
*         blank; data lines have the first digit of the MJD at the 
*         start.
          if( ierr.eq.0 .and. line(1:1).ne.' ' ) then
              read(line,120,iostat=jerr)  mjd, ixpr(1),ixpr(2),
     .                    iutr(1), iutr(2),ispr(1),ispr(2),
     .                    isut(1), isut(2),
     .                    num_sites,nf,num_svs,
     .                    ixpr(3),ixpr(4),
     .                    ispr(3),ispr(4), rho
 120          format(F8.2,2I8,I9,I7,2I6,2I8,I4,2I3,4I7,2x,3F7.3)
              call report_error('IOSTAT',jerr,'decod',line,1,
     .             'extend_iers/input file')
              num_ent = num_ent + 1
              if( num_ent.gt.max_ent ) then
                  write(*,*) 'Too many lines; Max allowed ',max_ent
                  stop 'Too many line in IERS file'
              end if

*             Now save the values that we need
              nm = num_ent
              if( nm.eq.1 ) then
*                 This is first entry; set the initial epoch
                  inf(1) = mjd
              end if
              if( nm.eq.2 ) then
*                 Get the spacing of the table
                  inf(2) = mjd - inf(1)
              end if
*             Check that the spacing is uniform
              djd = mjd-(inf(1)+inf(2)*(nm-1))
              if( abs(mjd-(inf(1)+inf(2)*(nm-1))).gt.0.001d0 ) then
                  write(*,170) nm, mjd, inf(1), inf(2)
 170              format('Entry ',I4,' at MJD ',F8.2,' not uniform',
     .                   ' for start at MJD ',F8.2,' with spacing ',
     .                   F4.2,' days')
                  stop 'extend_iers: IERS file not uniformly spaced'
              end if

*             Now save values
 
              pts(1,nm) = ixpr(1)/1.d6
              pts(2,nm) = ixpr(2)/1.d6
*             Get the number of leap seconds
              jd = mjd + 2400000.5d0
              call get_leapsec(jd, 1.d5, tai_utc)
              write(*,*) jd, tai_utc 
              pts(3,nm) = iutr(1)/1.d7 + tai_utc
*
*             Save the rates in case we need to use these later.
*             Convert LOD to dUT/dT
              rts(1,nm) =  ixpr(3)/1.d6
              rts(2,nm) =  ixpr(4)/1.d6
              rts(3,nm) = -iutr(2)/1.d7 

           end if
      end do
      inf(3) = num_ent
      
      write(*,*) inf
      do i = 1, num_ent
         write(*,*) (pts(j,i), j = 1,3)
      end do 
        


****  We have now read all of the IERS table into memory
      if( num_ent.gt.4 ) then

*         First create the pmu.<root> file.  This is a simple
*         extrapolation of the tables
          write(200,210) in_iers(1:trimlen(in_iers))
 210      format('* GLOBK PMU file created from ',a,
     .           ' with Lagrange extrapolation')
          do i = -2, num_ent+2
             mjd = inf(1) + inf(2)*i
             call lagrange_intp( inf, pts, mjd, value, 3)
             call jd_to_ymdhms( mjd+1.d-5, date, sec_tag)
             write(200,220) date, (value(j,1), 0.01d0, j=1,3)
 220         format(1x,i4,4(1x,i2),1x,4f9.6,F12.6,F9.6)
          end do

****      OK, Now write the GAMIT tables 
          st_mjd = nint(inf(1)) - 3
          en_mjd = nint(inf(1)+inf(2)*(num_ent+3))
*         Set number of entries per line
          npr = 6
          ut1_mult = 1.d-6 
          ut1_format = '(5x,i5,6(i9,1x),14x,i2)'
          pol_mult = 1.d-5
          pol_format = '(5x,i5,12i6,8x,i2)'


          write(201,250)  in_iers(1:trimlen(in_iers)),
     .           ut1_format, 4, nint(st_mjd+2400000.d0), 
     .           nint(en_mjd+2400000.d0), 
     .           npr, nint(inf(2)), ut1_mult

 250      format('UT1 values from IERS table file ',a,/,
     .           a32,1x,i1,1x,i7,1X,I7,2X,I1,2X,I1,10x,e7.1)

          write(202,255)  in_iers(1:trimlen(in_iers)),
     .           pol_format, nint(st_mjd+2400000.d0), 
     .           nint(en_mjd+2400000.d0), 
     .           npr, nint(inf(2)), pol_mult

 255      format('Pole values from IERS table file ',a,/,
     .           a35,i7,1X,I7,1X,I2,1X,I2,9x,e7.1)

C         do bjd = st_mjd, en_mjd, npr*inf(2)
          do ibjd = 0, nint((en_mjd-st_mjd)/(npr*inf(2)))
              bjd = st_mjd + ibjd*npr*inf(2)
              do i = 0,npr-1
                 mjd = bjd + inf(2)*i
                 call lagrange_intp( inf, pts, mjd, value, 3)
                 igxp(i+1) =  value(1,1)/pol_mult
                 igyp(i+1) =  value(2,1)/pol_mult
                 igut(i+1) = -value(3,1)/ut1_mult
              end do
              write(201,ut1_format) nint(bjd), (igut(i),i=1,npr)
              write(202,pol_format) nint(bjd), 
     .                              (igxp(i), igyp(i), i=1,npr)
              write(*,*) nint(bjd), (igut(i),i=1,npr)
              write(*,*) nint(bjd), 
     .                              (igxp(i), igyp(i), i=1,npr)
           end do

****       Thats all for the long extrapolation
           close(200)
           close(201)
           close(202)
      end if

      end



  



