      program decyr_atm
** written by Kurt Feigl
** mods for day of year option M.Burc Oral
** mods for day of year option M.Burc Oral  Fri Jun  3 17:49:09 EDT 1994

* needs GAMIT idoy, julday

c     given SITE yy mm dd hh mm val sig
c     return SITE yy.yyyyy val sig

      integer iyr,imo,idy,ihr,imn,idoy
      real*4  val,sig, dyear

      integer*4 n_argc,iargc,ios
      character*4 type,site

      n_argc = iargc()

* read input filename from command line input
      call rcpar(1, type)


*** n_argc = 0   type is as defined by   kurt   year.decimalday
      if ( n_argc .eq. 0 )  type = "year"

*orig
*par   type site site#YY MM DD HH MM  ATM         UNC
*ATM_ZEN R DS1B  3    94  3  2 19 12  0.0806 +-   0.0042
*filtered
*DS1B     94  3  2 19 12  0.0806   0.0042

      do while (.true.)

         read (unit = 5,
     .      fmt     = *,
     .      iostat  = ios,
     .      err     = 1010,
     .      end     = 1020) site,iyr,imo,idy,ihr,imn,val,sig

           dyear =  idoy(iyr,imo,idy)  + (ihr + imn/60. ) /24.
c         write(*,*) site,iyr,imo,idy,ihr,imn,val,sig

c     deal with leap years
         if (mod(iyr,4) .eq. 0) then
           if(type.eq."year") dyear= iyr +  dyear /366.
         else
           if(type.eq."year") dyear= iyr +  dyear /365.
         endif

         write(*,'(1x,A4,3x,f10.3,2f15.4)') site,dyear,val,sig
      enddo


 1010 continue
      call ferror (ios,6)

 1020 continue

      end




