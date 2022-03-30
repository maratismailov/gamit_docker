      subroutine read_glbk_dat(ifil1,ifil2)
c
c     read GLOBK data file and create FONDA format data file
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name1,sitnam(120)
      character*30 note
      character*1 symb
      integer iyear,id,id1
      integer ifil1,ifil2,i,jobs,itp,igp,j 
      integer nblen,lcmd             
      integer*4 julday
      character*80 wcmd
c
c     default format
      if (infmt(1:1).eq.'*')
     . infmt(1:41) = '(a1,2(f8.3,1x),6f8.2,f7.3,2x,3f8.2,1x,a8)'
c     
c     First, sort out all working sites and experiment number.
      nsit = 0
      cor3 = 0.0001d0

      call getcmd(ifil1,'Solution refers',wcmd,lcmd,1)
      if (lcmd.gt.0) then
         j = nblen(wcmd)
         note(1:12) = wcmd(j-31:j-23)
         read (note,'(f9.4)') time1
         print*,' get reference time:',time1
      else
         print*,' Can not get reference time.'
      endif

      do 20 i = 1,10000
         read (ifil1,fmt=infmt,end=50,err=20) 
     .      symb,a1,a2,v1,v2,a3,a4,s1,s2,cor,v3,a5,s3,name1 
c        skip comment line
         if (symb.eq.'*') goto 20
            nsit = nsit+1
            sitnam(nsit) = name1
 20   continue
 50   print*,'nsit =',nsit
      if (nsit.lt.1) goto 120
      jobs = nsit
c
c     site id
      note = '{network site number          '
      write (ifil2,30) nsit,note
      write (ifil2,'(a8)') (sitnam(i),i=1,nsit)
c
c     experiment number (always 1)
      i = 1
      itp = 24
      note = '{experiment number            '
      write (ifil2,30) i,note
      note = '{exp. index, obs. type, number'
      write (ifil2,140) i,itp,jobs,note

c     second, shift the format one by one
      rewind (ifil1)
      igp = 1
      do 90 i = 1,10000
 85      read (ifil1,fmt=infmt,end=120,err=90) 
     .      symb,a1,a2,v1,v2,a3,a4,s1,s2,cor,v3,a5,s3,name1 
c        skip unrealistic data
c        skip comment line
         if (symb.eq.'*') goto 85
c        identify site
         do 70 j = 1,nsit
            if (name1.eq.sitnam(j)) then
               id1 = j
               goto 100
            endif
 70      continue

 100     write (ifil2,150) time1,v1,s1,v2,s2,v3,s3,cor,cor3,cor3,name1

 90   continue

c
c10   format (a1,f7.3,1x,f7.3,6f8.1,f7.3,3f8.1,2x,a8)
c10   format (a1,f7.3,1x,f7.3,6f7.1,f8.3,3f7.1,1x,a8)
c10   format (a1,f8.3,1x,f8.3,1x,6f7.1,f7.3,2x,3f7.1,1x,a8)
 30   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 150  format (f9.4,3(1x,f11.4,1x,f8.3),3f9.5,2(1x,a8))
c
 120  continue
      return
      end

