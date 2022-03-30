CTITLE
 
      subroutine plts2
c
c     Second part of the plot program.. This one does all
c     of the plotting
c
c Include files
c -------------
*                         ! the PLOT parameter file
      include 'plot_param.h'
c
*                         ! the PLOT_COM common block
      include 'plot_com.h'
c
*                         ! the ema declaration
      include 'plot_ema.h'
 
c Variable declarations
c ---------------------
c
c Functions used
c --------------
c Trimlen -- HP utility
c rcpar   -- HP runstring utility
c ierr    -- IOSTAT error
c iel     -- Program control
 
      integer*4 trimlen, ierr, iel
 
c
c     Start reading from the userlu
c     =============================
c
      iel = 0
c
*                             ! loop until iel is set to one
      do while (iel.ne.1)
c                               which will occurr when the end
c                               command is issued
c....    Get next command
*                                   ! clear the buffer
         buffer = ' '
c
c....    Send a prompt if interactive
         if( terminal .and. pushlu.eq.0 ) then
            write(userlu,'(''?'',$)')
            read(*,'(a)', iostat=ierr) buffer
         elseif ( terminal ) then
            write(*,'(''KEY?'',$)')
            read(*,'(a)', iostat=ierr) buffer
         else
            read(userlu,'(a)', iostat=ierr) buffer
         end if
c
c....    See if EOF
         if( ierr.ne.0 ) then
              pcontrol = 2
              pel = 1
              iel = 1
              buffer = ' '
         end if
 
c....    If there is anything in the buffer and the first character is
c        blank (indicating a non-comment line) process it
         if( trimlen(buffer).gt.0 .and. buffer(1:1).eq.' ' ) then
c
c.....      Get the command
            call get_cmd(buffer, commands, num_commands, iel)
c
c....       Now process the command if one was found
*                                     ! read rest of command and process
            if( iel.gt.0 ) then
               pcontrol = 2
               call process( iel, ema_data(ibak_array),
     .            ema_data(ix_array),ema_data(iy_array),
     .            ema_data(ipt_array) )
c
*                                         ! this operation requires plot1
               if( pcontrol.ne.2  ) then
c
c....             Set iel to one to force exit after saving the command
c                 the program control version of iel
                  pel = iel
                  iel = 1
               end if
c
*                                     ! could not find command
            else
               call command_error(buffer,iel)
            end if
c
*                               ! echo comment to report device
         else
c
            call report(buffer)
c
         end if
c
      end do
c
c.... Thats all
c
*                   ! flush the buffer
      call jmcur
      return
c
      end
 
