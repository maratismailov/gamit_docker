      Subroutine GETSNAME( siteid,sname,sitnam,itype )

c    Move the characters of the site name to the right place for apr-style 
c    or l-file style output
                    
      integer*4 itype

c       itype = -2 : convert an 8-character apr-file site name to a 16-character site name for TFORM
c             = -1 : convert a 4-charcter L-file id and 12-character name into a 16-character 
c                      site name for TFORM
c             =  1 : convert a 16-character TFORM site name into a 4-character L-file id 
c             =  2 : convert a 16-character TFORM site name into an 8-character apr-file id

c       If the first character of any name is blank, we assume the input was via the 
c       terminal with no name; in this case keep it blank.
                 
      character* 8 siteid      
      character*12 sname
      character*16 sitnam



      if( itype.lt.0 ) then  

        sitnam = ' ' 
        if( siteid(1:1).ne.' ') then 
          sitnam(1:4)  = siteid(1:4)
          if( itype.eq.-2 ) then
            sitnam(5:12) = siteid(1:8)   
          else
            sitnam(5:16) = sname  
          endif
        endif
               
      elseif( itype.gt.0 ) then 
        siteid = ' ' 
        siteid(1:4) = sitnam(1:4)
        sname = ' '
        if( sitnam(1:1) .ne.' ') then   
          if( itype.eq.2 ) then
            siteid(5:8) = '_GPS'
          else
            sname=sitnam(5:16)
          endif
        endif

      endif

      return
      end 
               
   
