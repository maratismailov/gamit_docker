      program edit_apr
 
      implicit none 

*     This program will edit the apr file to round the coordinates to the
*     nearest number of meters specified in the input;  the purpose is to
*     allow distribution of apr files in cases where precise coordinates
*     of restricted by government policy. The runstring for the program is:
*
*     % edit_apr [edit file] [trunc]
*
*     where [edit_file] is an apr file containing
*     and   [roundoff] is a value in meters to which the coordinates are
*                      are to be truncated is a list of file to be operated on.  
*     The output of this program is currently to stdout, which can be redirected
*     to a file.
         
    
      character*8   site,ctrun
      character*16  edit_file
      character*256 line
      
      integer*4 ierr,len_run,trimlen,i,rcpar,itrun

      real*8 x(3),v(3),vsig(3),epoch,xtemp
 
****  Decode the runstring and open the files

      len_run = rcpar(1, edit_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'edit_apr', 'edit_apr', 1)
      end if               
      open(1, file=edit_file, iostat=ierr, status='old')
      if( ierr.ne.0 ) then
        call report_error('IOSTAT',ierr,'open',edit_file,1,
     .                'edit_apr')
      endif
      len_run = rcpar(2, ctrun)
      read(ctrun,*) itrun
      if( len_run.le.0 ) then
          call proper_runstring( 'edit_apr','edit_apr',1)  
      endif
      write(*,*) "Truncating coordinates to ",itrun," meters"
      if( itrun.eq.0 ) then
        write(*,*) "Input truncation = 0, assuming 1 m"
        itrun = 1
      endif      
cc      ii = trimlen(edit_file) 
cc      out_file = edit_file
cc      out_file(ii+1:ii+6) = ".trunc"                                              
cc      write(*,*) "Output file is ",out_file
cc      open(2,file=out_file,iostat=ierr,status='unknown')
cc      call report_error('IOSTAT',ierr,'open',out_file,1,
cc     .                  'out_file')
 
***** Now read the edit file
 
      do while ( ierr.eq.0 )
        read(1,'(a)', iostat=ierr ) line  
        if( ierr.eq.0 .and. trimlen(line).gt.0 ) then
          if( line(1:1).ne.' ' ) then
           write(*,'(a)') line
          else
           read(line(2:9),'(a8)') site
           read(line(10:256),*) x,v,epoch,vsig
           do i=1,3                
             xtemp = anint(x(i)/itrun) 
             x(i) = xtemp*itrun 
           enddo
           write(*,'(1x,a8,3(1x,f9.0),3(1x,f7.4),1x,f9.4,3(1x,f7.4))')
     .            site,x,v,epoch,vsig
          endif
        end if
      end do
 
****  Thats all
      end
 
 
 
 
 
 
 
 
 
 
