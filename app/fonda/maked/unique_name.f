      integer function unique_name(nsit,netid,sitnam,
     .                 abort_name,net_abt,id_list)
c
c     FONDA use 8-character site id to identify site name.
c     Unfortuntely, some site names of USGS data are identical 
c     even to 8-characters.  Unique-name technique must be exerted
c     to ensure unique 8-character site id.
c

      character*5 netid(nsit),nettmp(10),netido,net_abt(50)
      character*8 name1,sitnam(nsit)
      character*8 abort_name(50),namtmp(10)
      integer id_list(50)
      logical same
      integer i,j,isame,idup,inet,nsit
c     
c     First, sort out all working sites and find duplicated 
c     site names
      idup = 0
      do 20 i = 1,nsit-1
         name1 = sitnam(i)
         netido = netid(i)
         isame = 0
         same = .false.
         do j = i+1,nsit
            if (name1.eq.sitnam(j)) then
               same = .true.
               isame = isame+1
               nettmp(isame) = netid(j)
               namtmp(isame) = sitnam(j)
            endif
         enddo
c        find duplicated site name
         if (same) then
            idup = idup+1
            net_abt(idup) = netido
            id_list(idup)  = i
            abort_name(idup) = name1
            if (isame.eq.1) then
               sitnam(i)(6:8) = '_' // netido(1:2)
            else
               inet = 0
               do j = 1,isame
                  if (netido.eq.nettmp(j)) inet = inet+1
               enddo
               if (inet.lt.1) then
                  sitnam(i)(6:8) = '_' // netido(1:2)
               else
                  if (inet.eq.1) 
     .               sitnam(i)(6:8) = '_' // netido(1:1) // '0'
                  if (inet.eq.2) 
     .               sitnam(i)(6:8) = '_' // netido(1:1) // '1'
                  if (inet.gt.2)
     .               sitnam(i)(6:8) = '_' // netido(1:1) // '2'
               endif                 
            endif
            print*,' --- site name (',netido,') ',name1,
     .             ' has been changed to ',sitnam(i)
            goto 20
         endif
 20   continue
c
c     Theoretically, this loop is not enough.  Duplicated site name
c     still possibly exists.  For the time being, let's see if we
c     can get unique name list.
                  
      unique_name = idup            

      return
      end
