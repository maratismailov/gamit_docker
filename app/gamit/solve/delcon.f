      Subroutine DELCON( sitpos,numsit,outcon )
c
c     J.L. Davis 870308
c     mod:  J.L. Davis, M.H. Murray 870729    corrected bug
c           due to including all pairs of differences.
c           Changed so that only independent differences are used.
c
c     Routine to calculate the atmospheric-delay constraint matrix.
c     This is a matrix for constraining the atmospheric zenith-delay
c     parameter differences.
c
c     Input variables
c     ---------------
c     sitpos      A vector of geocentric cartesian site coordinates,
c                 in km (see module ATMCON).  The coordinates are
c                 stored as triplets of numbers, with the first in
c                 the triplet being the X-coordinate, the second Y,
c                 and the third Z.  For example, SITPOS(4) is the
c                 X-coordinate for the second site.
c
c     numsit      The number of sites in the solution.  It is assumed
c                 that there is one zenith delay per site.    
c

      implicit none

      include '../includes/dimpar.h'

      integer maxdif, dimop, condim, numsit
c
      parameter (maxdif = maxsit * (maxsit - 1) / 2)
      parameter (dimop  = 2 * maxdif)
      PARAMETER (CONDIM = MAXSIT * (MAXSIT + 1) / 2)

      integer newnum,i

      real*8 sitpos(MAXCRD)

      integer iusest,iusesa,iatcon
      common/stwght/iusest(maxsit),iusesa(maxsat),iatcon
c
c     Output variables
c     ----------------
c     outcon      This is the matrix which may be added directly to the
c                 normal matrix (i.e. A G A(T), where G is the weight
c                 matrix) after expansion.  The weights for the zenith
c                 delays are in cm^(-2).  OUTSIT is stored in lower
c                 triangular form, so its size is NUMSIT*(NUMSIT+1)/2.
c
      real*8 outcon(condim)
c
c     Parameters
c     ----------
c     maxsit      The maximum number of sites allowed.
c
c     maxdif      The maximum number of rows for the difference
c                 operator.  For MAXSIT sites, there are
c                 MAXSIT*(MAXSIT-1)/2 differences.
c
c     condim      The dimension of the original constraint matrix, in
c                 lower diagonal form.
c
c     dimop       The dimension for the "packed" operator.  Normally,
c                 the difference operator is a MAXDIF x MAXSIT 2-D
c                 array.  However, using the "Bock packing technique,"
c                 only the non-zero elements are stored in a 1-D
c                 array.  There are two non-zero elements for every
c                 difference (+1 & -1), so DIMOP will be 2 * MAXDIF.
c
c
c
c
c     Internal variables
c     ------------------
c     atmmat      The original matrix with the structure-function
c                 entries (see ATMCON)
c     difop       The difference operator
c     difopt      Transpose of difference operator
c     err         Error flag on inverting "temp" (not checked)
c     iwork       Work matrix
c     nrdif       Nr rows, difference matrix
c     ncdif       Nr columns, difference matrix
c     nrdift      Nr rows, transpose of difference operator
c     ncdift      Nr columns, transpose of difference operator
c     nrtemp      Nr rows, "temp" matrix (see below)
c     nrtemp      Nr columns, "temp" matrix
c     pntdf       The "pointer" array for the packed difference operator
c     rowdf       The "column counter" for the packed diff. operator
c     pntdft      "Pointer" array, transpose of diff. op.
c     rowdft      "Column pointer," transpose of diff. op.
c     temp        Temporary matrix for storage of DT A D (see below)
c     work        Working array
c
      integer
     .        err
     .    ,   iwork(maxsit)
     .    ,   nrdif
     .    ,   ncdif
     .    ,   nrdift
     .    ,   ncdift
     .    ,   nrtemp
     .    ,   pntdf(dimop)
     .    ,   rowdf(maxdif)
     .    ,   pntdft(maxdif)
     .    ,   rowdft(dimop)
c
      real*8
     .        atmmat(condim)
     .    ,   difop(dimop)
     .    ,   difopt(dimop)
     .    ,   temp(maxdif)
     .    ,   work(maxdif)
     .    ,   rcond
c
c     Equaivalences
c     -------------
c     Equivalence work arrays together to save space
c
      equivalence (work, iwork)
c
c     External references
c     -------------------
c     atmcon          Builds original atmospheric constraint matrix
c     atmdif          Forms the difference operator
c     atpa            Calculates matrix outer product
c     invers          Inverts square symmetric matrix
c     trans           Calculates the transpose of a matrix
c
c.... Form the matrix with original constraint matrix.
C---- Fit multi-session mode.   -DND- 880426
      call atmcon( sitpos,numsit,atmmat )
c
c.... What is the dimension of this square matrix?
      NEWNUM=0
      DO 100 I=1,NUMSIT
        IF (IUSEST(I).EQ.0) GO TO 100
        NEWNUM=NEWNUM+1
 100  CONTINUE
c
c.... Form the difference operator
      call atmdif(NEWNUM,difop,pntdf,rowdf)
c
c.... What are the number of rows and column for the difference matrix?
c     mod:  870729 see above
      nrdif = NEWNUM - 1
      ncdif = NEWNUM
c
c.... Form the product D A DT, where D is the difference operator, and
c     A is the constraint matrix.  We do this by using the routine
c     ATPA.  This routine, as its name implies, calculates AT P A,
c     but if we pass the transpose of the matrix which it needs, it
c     will calculate A P AT.  Since this subroutine expects the
c     transpose of A to be passed, we will pass the original matrix.
c
      call atpa(difop,pntdf,rowdf,nrdif,ncdif,atmmat,temp,work)
c
c.... What are the # of rows and columns of this "temporary" matrix?
      nrtemp = nrdif
c
c.... Invert this "temporary" matrix
      call inver2(temp,work,1,nrtemp,rcond,err)
      if (err.eq.130)
     .  call report_stat('WARNING','SOLVE','delcon',' '
     .     , 'Bad inverse of atmospheric constraint matrix',0)
c
c.... Determine the transpose of the difference matrix
      call trans(difop,difopt,pntdf,pntdft,rowdf,rowdft,nrdif,ncdif,
     .                                                             work)
c
c.... What are the dimensions of the transposed difference matrix?
      nrdift = ncdif
      ncdift = nrdif
c
c.... Now form the product DT (temp)^(-1) D
      call atpa(difopt,pntdft,rowdft,nrdift,ncdift,temp,outcon,work)
c
      end
