      Subroutine sort_string( ndim, in_num, in_names, new_num, new_names 
     .                      , out_num, out_names )
                        
c     Reorder a character array to add new entries, alphabetize, and 
c     remove duplicates.   Generalized from D. Dong subroutine ordsit,
c     called by CTOX, FIXDRV, XTORX, and now MAKEXP.   R. King 970307

         
c     Input:   ndim             I*4    Dimension of name arrays 
c              in_num           I*4    Number of non-blank elements in input name array
c              in_names(ndim)   CHAR   Input name array
c              new_num          I*4    Number of elements to be added to name array
c              new_names (ndim) CHAR   Array of new elements to be added to in_names
c              
c     Output   out_num          I*4    Number of non-blank elements in output name array
c              out_names(ndim)  CHAR   Output name array
               

      implicit none
    
      integer*4 ndim, in_num, new_num, out_num, live_num, mlen, slen
     .        ,  rcpar, i, j

      character*(*) in_names(ndim), new_names(ndim), out_names(ndim)
      character*256 prog_name, message, buf1, buf2

c     Get calling module name for report_stat 
      mlen = rcpar(0,prog_name)
               
c     Get length of strings
      slen = len(in_names(1))

c     First copy the input array to avoid overwriting it
                        
      do i = 1,ndim
        out_names(i) = ' '
      enddo
      do i = 1,in_num
        out_names(i) = in_names(i)
      enddo
      live_num = in_num
        

c     If new elements are to be added, do this before sorting
        
      if( new_num.gt.0 ) then
        do i = 1, new_num
           out_names(in_num+i) = new_names(i)
        enddo
        live_num = in_num + new_num
        if( live_num.gt.ndim ) then
           write(message,'(a,i4,a,i4,a)')   
     .     'Augmented array size (',live_num,') exceeds dimensions ('
     .     ,ndim,')'
           call report_stat('FATAL',prog_name,'makexp/sort_string'
     .                     ,' ',message,0)
        endif  
      endif
           

c     Now sort alphabetically

      do i = 1,live_num-1 
        do j = 1,live_num-i
           buf1(1:slen) = out_names(j)
           buf2(1:slen) = out_names(j+1) 
           if( lle(buf1,buf2) ) then
             out_names(j) = buf1(1:slen)
             out_names(j+1) = buf2(1:slen)
           else
             out_names(j) = buf2(1:slen)
             out_names(j+1) = buf1(1:slen)
           endif 
         enddo
      enddo

c     Finally, remove any duplicates
        
      i= 1
      do while (i.lt.live_num)  
         if ( out_names(i+1).eq.out_names(i) ) then
             do j = i,live_num-1
                out_names(j) = out_names(j+1)
             enddo 
             live_num = live_num - 1   
         else
             i = i + 1
         endif 
      enddo              
      out_num = live_num 
      do i = out_num+1,ndim
        out_names(i) = ' '
      enddo
      
      return
      end


