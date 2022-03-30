      subroutine getcmd(ifile,wkey,wcmd,icmd,job)
c
c     get command through key word approach
c
      character*256 line
      character*(*) wkey,wcmd
      integer i,icmd,ix,length,ifile,job,iy,mchkey,strjst
      integer nblen,ierr

c     ifile   : unit number for input file to be read

c     wkey    : character string of input command to be searched for

c     wcmd    : character string of command value, to be returned

c     icmd    : length of return string (wcmd)
c               = -1 input command not found
c               =  0 input command found but no value (ok in some cases)
c               >  0 value value returned      

c     job = 1 : rewind file and search for key word
c     job = 2 : rewind file and stop on key word
c     job = 3 : search for keyword without rewind.            
c     
      ierr = 0
c            
c     set wcmd length to -1 to indicate default wkey not found
      icmd = -1
      call blank(wcmd)
      call blank(line)
      iy = len(wkey)
c      print *,'wkey,job, ',wkey,job
      go to (50,100,30) job

c
c     search key word from the begining
 50   rewind(ifile)
 10   read(ifile,'(a256)',end=200,iostat=ierr) line  
      if (ierr .ne. 0 ) goto 200       
      if( line(1:1).ne.' ' ) goto 10
      ix = index(line,':')  
c     'end:' to terminate a command group used only in FONDA, not GAMIT; ok to
c     leave (I think) provided we properly skip comment lines (line added after 10 above)
c     and look for 'end:' rather then 'end' --rwk 980217
      i = mchkey(line,'end:',10,4)
      if (i.ge.0) goto 200
c     skip too short line 
      if (ix.le.iy) goto 10
c     skip comment line
      if (line(1:1).ne.' ') goto 10
c     skip mismatch line
      i = mchkey(line,wkey,ix,iy)
      if (i.lt.0) goto 10
      length = nblen(line)
c     blank command -- set length to zero
      if ( length.le.ix ) then
       icmd = 0
       goto 200  
      else
        icmd = length-ix
        wcmd(1:icmd) = line(ix+1:length)
        icmd = strjst(wcmd)
        icmd = nblen(wcmd)
        goto 200
      endif
c
c
c     search key word and stop there
 100  rewind(ifile)
 20   read(ifile,'(a256)',end=200,iostat=ierr) line
      if (ierr .ne. 0 ) goto 200      
      ix = index(line,':')
      iy = len(wkey)
      i = mchkey(line,wkey,ix,iy)
      if (i.lt.0) goto 20
      if (i.gt.0) then
        icmd = 0 
        goto 200
      endif
c
c     find the key word line without rewind
 30   read(ifile,'(a256)',end=200,iostat=ierr) line
      if (ierr .ne. 0 ) goto 200
      ix = index(line,':')
      iy = len(wkey)
c     exit     
c     print *,'line wkey ix iy ',line(1:30),wkey,ix,iy,ierr
      i = mchkey(line,'exit',10,4)
      if (i.ge.0) goto 200
c     skip too short line
      if (ix.le.iy) goto 30
c     skip comment line
      if (line(1:1).ne.' ') goto 30
c     skip mismatch line
      i = mchkey(line,wkey,ix,iy)
c      print *,'wkey ix iy i ',wkey,ix,iy,i
      if (i.lt.0) goto 30
      length = nblen(line)
c     blank command - set length to zero
c      print *,'ix length ',ix,length
      if ( length.le.ix ) then
        icmd = 0
        goto 200
      else
        icmd = length-ix
        wcmd(1:icmd) = line(ix+1:length)
        icmd = strjst(wcmd)
        icmd = nblen(wcmd) 
c        print *,'wcmd icmd ',wcmd,icmd
      endif
c
 200  continue
      return
      end
