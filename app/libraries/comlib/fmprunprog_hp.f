CTITLE FMPRUNPROGRAM
 
      subroutine fmprunprogram( name, ifive, runname, num_save, 
     .                          save_unit, offset_unit )

      implicit none
 
*     This is a sub routine which will simply tell the user the
*   next program needs running.  In the UNIX version, running the
*     programs will probably be replaced with subroutine calls.
*
* variables
*   name        - Name of program to run with its
*               - runstring
*  runname      - Nomimals the actual run name
 
      character*(*) name, runname
 
*   ifive(5)    - Five element vector with is not used.
*   num_save    - Number of units to be saved
*   save_unit(num_save) - Unit numbers whose pointers are to be
*               - saved
*   offset_unit(num_save) - Current pointer offset from the begining
*               - of the file.  Value is restored when the system
*               - call is completed
 
      integer*4 ifive(5), num_save, save_unit(num_save), 
     .          offset_unit(num_save)
 
* Local varibales
 
*   trimlen     - Used for length of string
*   i - loop counter
*   iru - position of ru command
*   ierr   - Fseek error return
*   system - SUN4 function to shell a command
*   Ftell  - Returnes pointer position for file (Sun4 library)
*   Fseek  - Positions file to offset location (value obtained
*          - Ftell call)
 
      integer*4 trimlen, i, iru, ftell, ierr
* PT & HM 960920: make "system" an external variable for g77 compilation
* TAH 061220: Needed for HP version 
      external system
      integer*4 system

***** Write out the program to run
*     Remove any commas from the names
      do i = 1, trimlen(name)
          if( name(i:i).eq.',' ) then
              name(i:i) = ' '
          end if
      end do

***** Remove the RU from the front
      iru = index(name, 'RU')
      if( iru.gt.0 ) then
          name(iru:iru+1) = ' '
      end if

      write(*,100) name(1:max(1,trimlen(name)))
 100  format(/,' Program ',a,' running ')

****  Save the offsets of the files which are desired to be 
*     kept open

      do i = 1, num_save
         offset_unit(i) = ftell(save_unit(i))
      end do

***** NOW shell out to command
      ifive(1) = system(name)
 
****  If one retunred then stop
      if( ifive(1).ne.0 ) then
          stop ' Program stopped in FmpRunProgram'
      end if

***** Now restore the offsets to the file
      do i = 1, num_save
         if( offset_unit(i).lt.0 ) then
             write(0,200) save_unit(i), offset_unit(i)
  200        format(' **** FTELL ERROR **** on file unit ',i5,
     .              ' Ftell returned offset was ',i10,/,
     .              ' Error reported by FmpRunProgram')
         else
*                                                 ! 0 says position from
*                                                 ! start of file.
c** function     ierr = fseek(save_unit(i), offset_unit(i), 0)
c** intrinsic    call fseek(save_unit(i), offset_unit(i), 0)
             call fseekg(save_unit(i), offset_unit(i), 0, ierr)
             call report_error('PERROR',ierr,'fseek',
     .           'file unit after SYSTEM call',0,'FmpRunProgram')
         end if
      end do
 
****  Thats all
      return
      end
 
