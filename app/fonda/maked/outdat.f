      subroutine outdat(idatf,imapf)
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*30 note
      integer idatf,imapf,i,is,igo,jobs,ie,itp,igp
      integer it1,it2,ig,iobs,ib,kobs,id1,id2,ide,imi
      integer*4 julday
c
c     output mapping file
      if (iomode(5).gt.0) then
c     header:
      write (imapf,'(2x,''velocity field'')')
c     data
      do 20 i = 1,nsit
         is = itoj(i)
         if (comode.le.2) then
            alon = slo(is)*rtod
            alat = sla(is)*rtod
            write (imapf,10) alon,alat,ve(is),vn(is),sname(is)
         else
            write (imapf,10) x(is),y(is),vx(is),vy(is),sname(is)
         endif
 20   continue
      endif
 10   format (1x,2f11.5,2f9.3,2x,a8)
c
c     output data file
      igo = 0
      jobs = 0
c     temporary
      rho = 0.0d0
c
c     site id
      note = '{network site number          '
      write (idatf,30) nsit,note
      write (idatf,'(a8)') (sname(itoj(i)),i=1,nsit)
c
c     experiment number
      note = '{experiment number            '
      write (idatf,30) iexp,note
 30   format (i5,25x,a30)
      do 50 ie = 1,iexp
         itp = itype(ie)
         note = '{exp. index, obs. type, number'
         write (idatf,40) ie,itp,inobs(ie),note
c        group number and time for each experiment
         igp = igroup(ie)
         it1 = itime(ie,1)
         it2 = itime(ie,2)
         time1 = 1900.0d0+julday(1,it2,it1,2)/365.2422d0
         do 60 ig = 1,igp
            iobs = igobs(ig+igo)
            do 70 ib = 1,iobs
            kobs = jobs+ib
c           observable and index
            id1 = jtoi(idobs(kobs,1))
            id2 = jtoi(idobs(kobs,2))
            if (itp.le.3) then
            call rtodms(data(kobs),ide,imi,s1,2)
            write (idatf,120) 
     *      time1,ide,imi,s1,erd(kobs),sname(itoj(id1)),sname(itoj(id2))
            elseif (itp.gt.3.and.itp.le.10) then
            write (idatf,80) 
     *      time1,data(kobs),erd(kobs),sname(itoj(id1)),sname(itoj(id2))
            elseif (itp.gt.20.and.itp.le.24) then
            write (idatf,110) 
     *      time1,data(kobs),erd(kobs),rho,sname(itoj(id1))
            elseif (itp.gt.24.and.itp.le.30) then
            write (idatf,110) 
     *      time1,data(kobs),erd(kobs),rho,sname(itoj(id1))
            else
            write (idatf,110) 
     *      time1,data(kobs),erd(kobs),rho,sname(itoj(id1))
            endif
 70         continue
            jobs = jobs+iobs
 60      continue
         igo = igo+igp
 50   continue
 40   format (1x,3i5,14x,a30)
 80   format (1x,f10.4,f14.5,f14.8,2x,a,2x,a)
110   format (1x,f10.4,f14.5,f14.8,f9.4,2x,a,2x,a)
120   format (1x,f10.4,2i5,f14.8,f9.4,2x,a,2x,a)
c
      return
      end

