      integer function mchkey(strinx,striny,ix,iy)

c     Find a string within another string
c
c     Input:
c        stringx:  string to be searched
c        stringy:  candidate string to be found
c        ix     :  length of stringx    
c        iy     :  length of stringy
c           
c     Output:  Value of mchkey (-1 if not found, +1 if found)
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

