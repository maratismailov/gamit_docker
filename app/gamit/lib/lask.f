      logical function lask( )
c     ask user yes or no question
      implicit none
      character*1 yesno,lowerc
      character*80 buff80, prog_name
      integer nerr,len,rcpar
c

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)

      nerr = 0 
      lask = .false.
  50  continue

      if (nerr .lt. 5) then
         write (6,10)
  10     format (1x,'(Y/N) ',$)
         read (5,20) buff80
  20     format (a)

         call ljust (80,buff80)
         yesno = buff80(1:1)

         if (lowerc(yesno) .eq. 'y') then
            lask = .true.
         else if (lowerc(yesno) .eq. 'n') then
            lask = .false.
         else
            write (6,30)
 30         format (/,1x,'Inappropriate response.',$)
            nerr = nerr + 1
            goto 50
         endif
      else
         call report_stat('FATAL',prog_name,'lib/lask',' '
     .                      ,'Tried 5 times',0)
      endif

      return
      end
