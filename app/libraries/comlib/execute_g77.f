ctitle
 
      subroutine execute( prog_name, ifive, num_save, save_unit,
     .                    offset_unit )
 
      implicit none

c     routine to exucute a program immediately with a wait by the
c     scheduling program.
c
c MOD TAH 870108 Modifed program to use FmpRunProgram instead of
c     XQPRG.  The new version should run either CI or FMGR programs
c
c Variables
c ---------
c prog_name -- the name of the program to be scheduled
c ifive -- option five paremeter rmpar thish can be passed to the
c     scheduled program.  These vals will be replaced by the
c     five paramters passed back from the scheduled program
* num_save  - Number of units to save positions of after call to
*           - system
* save_unit - Unit numbers to save
* offset_unit - Actuall offsets of the save units (see Fseek, and Ftell in
*               Sun4 Fortran manual.)
c
 
      integer*4 ifive(5), num_save, save_unit(num_save),
     .          offset_unit(num_save)
 
      character*(*) prog_name
 
c
c Local Variables
c ---------------
 
*   ierr            - Error returned from FmpRunProgram
*   FmpRunProgram   - HP subroutine (used to be a function)
 
      integer*4 ierr
 
*   RunName         - Actual name of the running program, i.e,
*                   - users LU number is added.
c**      character*24 RunName

c      integer*4 system  ! System routine to span process
c      external system 
 
c
***** Run the program
 
c**      call FmpRunProgram( prog_name, ifive, RunName, num_save,
c**     .                      save_unit, offset_unit  )

***** NOW shell out to command
      ifive(1) = system(prog_name)

      ierr = ifive(1)
 
c.... see if an error ocurred
*                           ! error did ocurr
      if( ierr.ne.0 ) then
 
         call report_error('FmpRunProgram', ierr,'run',prog_name,0,
     .                     'EXECUTE')
 
      end if
c
c.... That's all
      return
      end
 
