CTITLE add_GGV

      subroutine add_GGV( gdescription )

      implicit none

*     Routine to add last line of GGVerion file to
*     solution description

* PASSED
      character*(*) gdescription

* LOCAL
      character*128 home   ! User home directory
      character*128 ggv_file   ! Name of file with GGVersion
      character*128 line   ! line read from GGV file

      integer*4 ierr   ! IOSTAT error
      integer*4 unit   ! Unit for reading GGV files
      logical unit_open   ! True if unit number open

****  Get the home directory and create GGV file name
      call getenv('HOME',home)
      ggv_file = trim(home) // '/gg/com/GGVersion' 

*     Get unit to use
      unit_open = .true.
      unit = 150
      do while ( unit_open .and. unit.lt.300 )
        unit = unit + 1
        inquire( unit=unit, opened=unit_open )
      enddo
      if( unit.eq.300 ) then
        write(*,120) unit
 120    format('** WARNING ** No open units between 150 and ',i4)
        RETURN
      endif

****  Now read to end of file
      open(unit,file=ggv_file,status='old',iostat=ierr)
      if( ierr.eq.0 ) then
         do while ( ierr.eq.0 ) 
            read(unit,'(a)',iostat=ierr) line
         enddo
* Format if line from GGVersion
!10.71.006 Fri Mar 20 08:56:51 EDT 2020 Release version, plus extend stinfo comment length.
*        Update gdescription.
         gdescription = trim(gdescription) // ' GGVer ' // line(1:38)
      endif

***** Thats all
      return
      end


        
         
        
