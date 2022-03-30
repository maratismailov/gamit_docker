      subroutine pickfns(list,file,nfil)

c     pick a file name from the available choices

      implicit none

c     UNIX wild card
      character*(*) list
      character*(*) file(10)
      character*80 files(1000),buff80
      character*1 letter
      character*20 form 
      character*256 message
      integer nch,len,i,nfil,ib
      integer ipick,nfiles
      integer ioerr,nerr,count_arg,lift_arg
      logical fcheck

      nerr = 0
      call getdir (list,1000,files,nfiles)

      if (nfiles .gt. 0) then
         write (6,'(/,a)') 'Available files:'
         do  i = 1,nfiles
            write (*,'(i3,1x,a)') i,files(i)
         enddo

         write (6,'(/,1x,a)') 'Enter file names or pick numbers:'
         nch = len(file(1))
         read (5,'(a80)') buff80

         nfil = count_arg(buff80)

         do  i = 1,nfil
            ib = lift_arg(buff80,form,i)
            letter = form(1:1)
            if (lge(letter,'1') .and. lle(letter,'9')) then
              read (form(1:ib),*,iostat=ioerr) ipick
              if (ioerr .eq. 0) then
                if (ipick .lt. 1 .or. ipick .gt. 1000) then 
                  write(message,'(a,i4)') 'ipick overflow, ipick=',ipick
                  call report_stat('FATAL','MAKEJ','pickfns',' '
     .                            ,message,ioerr)
                else
                     file(i) = files(ipick)(1:nch)
                endif
              else
                  call report_stat('FATAL','MAKEJ','pickfns',' ',
     .            'Error, i/o error getting list ',ioerr)
              endif
            else
               file(i) = form(1:ib)
            endif
c           double check existence of file
            if (.not.fcheck(file(i))) then
              call report_stat('FATAL','MAKEJ','pickfns',file(i),
     .        'Error, opening file: ',ioerr)
            endif
         enddo
      else
         call report_stat('WARNING','MAKEJ','pickfns',list,
     .   'Cannot find: ',0)
      endif

      return
      end


