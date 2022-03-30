      integer function unique_nam1(nsit,iesit,sitnam)
c
c     FONDA use 8-character site id to identify site name.
c     Unfortuntely, some site names of USGS data are identical 
c     even to 8-characters.  Unique-name technique must be exerted
c     to ensure unique 8-character site id.
c

      character*8 name1,sitnam(nsit)
      character*8 abort_name(50),namtmp(50)
      integer id_list(50),iesit(nsit),nsit_abt(50)
      logical same
      integer i,j,j1,isame,idup,nsit,idsito
c     
c     First, sort out all working sites and find duplicated 
c     site names
      idup = 0
      do 20 i = 1,nsit-1
         name1 = sitnam(i)
         idsito = iesit(i)
         isame = 0
         same = .false.
         do j = i+1,nsit
            if (name1.eq.sitnam(j)) then
               same = .true.
               isame = isame+1
               write (namtmp(isame),'(i3)') iesit(j)
               if (iesit(j).lt.100) namtmp(isame)(1:1) = '0'
               if (iesit(j).lt.10) namtmp(isame)(2:2) = '0'
               id_list(isame) = iesit(j)
               nsit_abt(isame) = j
            endif
         enddo
c        find duplicated site name
         if (same) then
            do j = 1,isame
               idup = idup+1
               abort_name(idup) = name1
               j1 = nsit_abt(j)
               sitnam(j1)(6:8) = namtmp(j)(1:3)
               print*,' --- site name : ',name1,' (id = ',id_list(j),
     .            ')  has been changed to ',sitnam(j1)
            enddo
         endif
 20   continue
c
c     Theoretically, this loop is not enough.  Duplicated site name
c     still possibly exists.  For the time being, let's see if we
c     can get unique name list.
                  
      unique_nam1 = idup            

      return
      end
