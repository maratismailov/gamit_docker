      Subroutine ATMCON( sitpos,numsit,atmmat )
c
c     J.L. Davis 870305
c
c     Routine to set up atmospheric zenith-delay constraint matrix
c
c     Input variables
c     ---------------
c     sitpos          3-D site coordinates (units = km), arranged in
c                     a 1-D array with triplets (X-Y-Z) of coordinates.
c                     For example, SITPOS(1) is the X-coordinate for the
c                     first site, SITPOS(4) is the X-coordinate for
c                     the second site, and so on.
c
c     numsit          The number of stations.  The dimension of the
c                     atmospheric constraint matrix will be
c                     NUMSIT x NUMSIT.
c          
      implicit none

      include '../includes/dimpar.h'
c
      integer numsit, condim,jel,ij,i1,i2,j1,j2
c
c     Output variables
c     ----------------
c     atmmat          The atmospheric constraint matrix.  Since
c                     ATMMAT(I,J) = ATMMAT(J,I), the matrix is arranged
c                     in lower triangular form.  Thus, the element
c                     of ATMMAT with indices (I,J) is stored in
c                     ATMMAT(K), where K = I * (I - 1) / 2 + J, for
c                     I > J.
c
      PARAMETER (CONDIM = MAXSIT * (MAXSIT + 1) / 2)
      real*8 sitpos(MAXCRD)
      real*8 atmmat(condim)
c
c     Internal variables
c     ------------------
c     blen            Baseline length
c     i, j, ij        Indices
c     starti          Index in SITPOS of X-coord for ith site
c     startj          Index in SITPOS of X-coord for jth site
c
      real*8 blen
c
c.... in order to fit the station option, increase iusest to judge the
c.... live stations.   -dnd- 870826
      integer*4 IUSEST(MAXSIT),IUSESA(MAXSAT) 
      integer*4 iatcon
      COMMON/STWGHT/IUSEST,IUSESA,IATCON
c
      integer i, j, starti, startj
c
c     Subroutines or  functions called
c     --------------------------------
c     baslin          Calculates baseline length
c
c     atmstr          Calculates the atmospheric structure function
c
c     jel             Returns inex for lower triangular storage
c
      real*8 baslin, atmstr
c
c.... Loop over all sites
        I1=0
      do 100 i = 1, numsit
        IF (IUSEST(i).EQ.0) GO TO 100
        J1=I1
        I1=I1+1
C---- To fit multi-session mode.  -DND- 880426
c*         I2=I+ICFILE
c changed w/ removal of mulit-lession RWK 190506
         i2 = i 
c
c....     Calculate the index in SITPOS for the X-coordinate of this
c         site
          starti = 3 * (i2 - 1) + 1
c
c....     Loop over sites in rows below this
          do 200 j = i, numsit
        IF (IUSEST(j).EQ.0) GO TO 200
        J1=J1+1
         j2=j
c
c....         Calculate the index in SITPOS for the X-coord for this
c             site
              startj = 3 * (j2 - 1) + 1
c
c....         Calculate the length for baseline between sites I and J
              blen = baslin(sitpos(starti),sitpos(startj))
c
c....         Determine the index in ATMMAT
              ij = jel(i1,j1)
c
c....         Calculate the entry based on the structure function
              atmmat(ij) = - 0.5D0 * atmstr(blen)
c
  200     continue
c
  100 continue
c
      end
