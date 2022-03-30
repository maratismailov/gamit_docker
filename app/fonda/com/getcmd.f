      subroutine getcmd(ifile,wkey,wcmd,icmd,job)
c
c     get command through key word approach
c
      character*128 line
      character*(*) wkey,wcmd
      integer i,icmd,ix,length,ifile,job,iy,mchkey,strjst
      integer nblen

c     job = 1 : rewind file and search for key word
c     job = 2 : rewind file and stop on key word
c     job = 3 : search for keyword without rewind.

      icmd = -1
      call blank(wcmd)
      call blank(line)
      iy = len(wkey)
      go to (50,100,30) job
c
c     search key word from the begining
 50   rewind(ifile)
 10   read(ifile,'(a128)',end=200) line
      ix = index(line,':')
      i = mchkey(line,'end',10,3)
      if (i.ge.0) goto 200
c     skip too short line
      if (ix.le.iy) goto 10
c     skip comment line
      if (line(1:1).ne.' ') goto 10
c     skip mismatch line
      i = mchkey(line,wkey,ix,iy)
      if (i.lt.0) goto 10
      length = nblen(line)
c     blank command
      if (ix.ge.length) goto 200
      icmd = length-ix
      wcmd(1:icmd) = line(ix+1:length)
      icmd = strjst(wcmd)
      icmd = nblen(wcmd)
      goto 200
c     
c
c     search key word and stop there
 100  rewind(ifile)
 20   read(ifile,'(a128)',end=200) line
      ix = index(line,':')
      iy = len(wkey)
      i = mchkey(line,wkey,ix,iy)
      if (i.lt.0) goto 20
      if (i.gt.0) goto 200
c
c     find the key word line without rewind
 30   read(ifile,'(a128)',end=200) line
      ix = index(line,':')
      iy = len(wkey)
c     exit
      i = mchkey(line,'exit',10,4)
      if (i.ge.0) goto 200
c     skip too short line
      if (ix.le.iy) goto 30
c     skip comment line
      if (line(1:1).ne.' ') goto 30
c     skip mismatch line
      i = mchkey(line,wkey,ix,iy)
      if (i.lt.0) goto 30
      length = nblen(line)
c     blank command
      if (ix.ge.length) goto 200
      icmd = length-ix
      wcmd(1:icmd) = line(ix+1:length)
      icmd = strjst(wcmd)
      icmd = nblen(wcmd)
      goto 200
c
 200  continue
      return
      end
c
c------------------------------------------------------
      integer function mchkey(strinx,striny,ix,iy)
c
      character*(*) strinx,striny
      integer ix,iy,i,i1
c
      mchkey = -1
      i1 = 0
      do 10 i = iy,ix
         i1 = i1+1
         if (striny(1:iy).eq.strinx(i1:i)) goto 20
 10   continue
      return
c
 20   mchkey = 1
      return
      end
