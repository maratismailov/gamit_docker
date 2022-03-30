      Subroutine upnam3( file_in, file_out )

c     Update the file name for an apr file by either adding an 
c     underscore and a letter or, if an underscore is already present,
c     by changing the letter.  R. King 020807

      implicit none
       
      integer*4 icol  
      character*1 achar,lowerc
      character*(*) file_in, file_out

                                              
      file_out = file_in 
      icol = index(file_in,'.')
      if( file_in(icol-2:icol-2).eq.'_' ) then
c       file is already of the form _x.apr, just change the letter
        achar = file_in(icol-1:icol-1)
        call newchr(achar)
        file_out(icol-1:icol-1) = lowerc(achar)
      else
c       need to add an underscore and the letter 'a'
        file_out(icol:icol+5) = '_a.apr'
      endif

      return
      end
