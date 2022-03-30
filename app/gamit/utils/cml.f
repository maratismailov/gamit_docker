      program cml

c compare files line by line
c Kurt Feigl 5 March 87
c
      character*200 line1,line2
      character*80 file_name1,file_name2
      integer*4 ioerror1,ioerror2,icount,in1,in2


      icount = 1

      print *,'Start column (Max 200)'
      read *,in1
      if (in1 .le. 0 .or. in1 .gt. 200) in1 = 1

      print *,'Stop column (Max 200)'
      read *,in2
      if (in2 .le. 0 .or. in2 .gt. 200) in2 = 200

      print *, 'First file? '
      read '(a80)', file_name1
      print *, 'second file? '
      read '(a80)', file_name2

      open(unit = 10,
     +   file   = file_name1,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerror1)
      if (ioerror1.ne.0) then
         print *, 'Error opening file: ', file_name1
         goto 1000
      else
         print *, 'Opened file: ', file_name1
      endif

      open(unit = 11,
     +   file   = file_name2,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerror2)
      if (ioerror2.ne.0) then
          print *, 'Error opening file: ', file_name2
          goto 1000
      else
         print *, 'Opened file: ', file_name2
      endif

      print *, ' Found the following differences:  '
      print *, ' '

 10   continue
      read(unit   = 10,
     +   fmt    = '(a200)',
     +   err    = 1000,
     +   end    = 500,
     +   iostat = ioerror1)
     +   line1

      read(unit   = 11,
     +   fmt    = '(a200)',
     +   err    = 1000,
     +   end    = 500,
     +   iostat = ioerror2)
     +   line2

      if (line1(in1:in2) .ne. line2(in1:in2)) then
         write (6,110) icount
 110     format (1x,'Line',1x,i7)
         write (6,120) line1(in1:in2)
         write (6,120) line2(in1:in2)
 120     format (15x,a)
      endif

      icount = icount + 1
      if (ioerror1 .eq. 0 .or. ioerror2 .eq. 0) goto 10

 500  if (ioerror1.eq.-1) then
         print *, ' '
         print *, 'Reached end of file: ',file_name1
      endif
      if (ioerror2.eq.-1) then
         print *, 'Reached end of file: ',file_name2
      endif


 1000 if (ioerror1.ne.0 .and. ioerror1.ne.-1) then
         call ferror (ioerror1,6)
      else if (ioerror2.ne.0 .and. ioerror2.ne.-1) then
         call ferror (ioerror2,6)
      endif

      stop
      END

