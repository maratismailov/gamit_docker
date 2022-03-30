C
      Subroutine REMENU( nbias,ndex,last_nonbias )
C
C---- Reorder labels of independent bias parameters in Q-file menu.
C        Order :
C           L1, L2-L1 by length of baselines
C        -DND- 880106
C
C---- Parameter discription :
C        NBIAS -- Number of independent double-difference  L1 bias parameters.
C        IBAND --    IBAND=1     L1 mode
C                    IBAND=2     L1,L2 mode or LC mode
C        RLABEL -- Label of all parameters.
C        SITNAM -- Name of all stations.
C        ISPRN -- Index of all satellites.
C        NDEX -- Index of stations of every baseline   

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer maxbas
      parameter (maxbas=maxsit*(maxsit-1)/2)
              
      character*2 buf2
      character*4 sit4

      integer ndex(maxbas*2),nbias,last_nonbias,nb
     .      , ibod0,ib,i1,i2,i3,id10,id20,id1,id2,is1,is2
     .      , ipp,il,ii,jj,i

      real*8 dwv

      logical debug/.false./
                       
      if(debug) then 
        print *,'REMENU nbias last_nonbias ndex ',nbias,last_nonbias
        write(*,'(10i7)') (ndex(i),i=1,50) 
      endif    
      nb = nsite*nsat*iband
      ibod = 0
      ibod0 = ibod
      ibod = ibod + nbias*iband
c---- rewrite label for independent l1 bias
      ib = nbias
      do 40 i = 1,ib*iband
         i1 = i + ibod0
 40   half(i1) = 1

c---- loop through all biases
      if(debug) print *,'REMENU ib ',ib 
      do 10 i = 1,ib
         i3 = i + ibod0
         i1 = last_nonbias + i
         rlabel(i1)(1:1) = 'B' 
c*        this formerly held the session number
* MOD TAH 190615: Retain "1" as "session numner" for backwards compatabilty
*        with script that search biases and to avoid space in line that 
*        will affect awk column numbers.
         rlabel(i1)(2:2) = '1' 
c*         read(unit = buf1,fmt = '(a1)') rlabel(i1)(2:2)
         rlabel(i1)(3:5) = 'L1 '
         i2 = ipntdt(i)
         jj = (i2-1)*2
c---- determine two stations of the baseline
         id10 = ndex(jj + 1) 
         id20 = ndex(jj + 2)    
         sit4 =  cfiln(id10)(2:5) 
         call uppers(sit4)
         rlabel(i1)(6:9) = sit4
         rlabel(i1)(10:10) = '-' 
         sit4 = cfiln(id20)(2:5)
         call uppers(sit4)
         rlabel(i1)(11:14) = sit4
         rlabel(i1)(15:15) = ' ' 
         id1 = id10    
c         print *,'i1 jj ndex id10 cfiln rlabel '
c     .     ,i1,jj,ndex(jj),id10,cfiln(id10),rlabel(i1)
         jj = (i - 1)*4
c---- determine two satellites with the independent bias
         is1 = ipntd(jj + 1)
         id2 = ipnt2d(is1)
         is1 = id2 - (id1-1)*nsat
         write(unit = buf2,fmt = '(i2)') isprn(is1)
         read(unit = buf2,fmt = '(a2)') rlabel(i1)(16:17)
         rlabel(i1)(18:18) = '-'
         is2 = ipntd(jj + 3)
         id2 = ipnt2d(is2)
         is2 = id2 - (id1-1)*nsat
         write(unit=buf2,fmt='(i2)') isprn(is2)
         read(unit=buf2,fmt='(a2)') rlabel(i1)(19:20)

c---- calculate wide-lane biases
         if (iband.ge.2) then
         ipp = ibod0/2 + i
         ddwl(ipp) = wl0(id10,is1) - wl0(id10,is2)
     .           - wl0(id20,is1) + wl0(id20,is2)
         dwv =  vwl(id10,is1)**2 + vwl(id10,is2)**2
     .        + vwl(id20,is1)**2 + vwl(id20,is2)**2
         ddwv(ipp) = dsqrt(dwv)
         endif
c---- set wavelength factor for ambiguity resolution
c        narrowlane (L1 or L2)
         il = 1
         if( l2flag.eq.-2 ) il = 2
         if (lwave(id10,is1,il).eq.2.or.lwave(id10,is2,il).eq.2)
     *    half(i3) = 2
         if (lwave(id20,is1,il).eq.2.or.lwave(id20,is2,il).eq.2)
     *    half(i3) = 2
         if (iband.eq.2) then
         ii = i3 + nbias
         if (lwave(id10,is1,2).eq.2.or.lwave(id10,is2,2).eq.2)
     *    half(ii) = 2
         if (lwave(id20,is1,2).eq.2.or.lwave(id20,is2,2).eq.2)
     *    half(ii)= 2
         endif
c      print *,'REMENU iband,i3,ii,half ',iband,i3,ii,half
c      print *,'id10 id20 is1 lwave(2) ',id10,id20,is1
c     .       ,lwave(id10,is1,2),lwave(id20,is1,2)
 10   continue

c---- create labels for widelane biases
      if (iband.eq.1) goto 50
      ib = nbias*iband          
c      print *,'for WL nbias iband ib ',nbias,iband,ib
      do 20 i = nbias+1,ib
         i1 = last_nonbias + i
         i2 = i1 - nbias
         rlabel(i1)(1:2) = rlabel(i2)(1:2)
         rlabel(i1)(3:5) = 'L21'
         rlabel(i1)(6:20) = rlabel(i2)(6:20)       
c        print *,'i1 i2 rlabel ',i1,i2,rlabel(i1)
 20   continue
 50   i1 = last_nonbias + 1 + ib
      i2 = last_nonbias + nb  
c      print *,'Bias junk ib nb last_nonbias i1 i2 '
c     .       , ib,nb,last_nonbias,i1,i2
      do 30 i = i1,i2
         rlabel(i) = 'Bias junk           '  
c         print *,'i rlabel ',i,rlabel(i)
 30   continue
      return
      end
