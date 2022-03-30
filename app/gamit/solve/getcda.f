      Subroutine GETCDA(ifile,istat,ierrfl,omc,isnr,el)

c     get data record from c-file
c                            
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
                 
c      passed values
      integer*4 ifile,istat,ierrfl(maxsat)
      real*8 omc(maxdat,maxsat) 
                        
c      internal 
      integer*4 itepch,jsat,iprn,iyr,idoy,i,j,ii
      real*8  data(maxdat,maxsat),save(maxsav),el(maxsat),azim(maxsat)
     .     ,  sod,rclock,pres,temp,relhum
     .     ,  spare(maxspr),svcepc,svcl1  
      real*8  zendel,l1z,l1n,l1e,l2z,l2n,l2e,antaz
      real*4 ampl1(maxsat),ampl2(maxsat),atmlod4(3)
c       atmosphere and delay values stored as scalars since not used here
      real*8 elvdot,azmdot,nadang,atmdel,tau,drate,partial(2)
     .     , latr_sph,lonr,radius
      integer*4 mprn(maxsat),data_flag(maxsat)
     .        , isnr(maxdat,maxsat),okmet,nspare,nsave,msat,ndat
      character*4 sitecd

      logical debug/.false./
  
      if(debug) print *,'GETCDA 1 eopch a ',iepoch,(a(ii),ii=1,3)

c        read data records    
cd         if(iepoch.eq.953) print *,'GETCDA iepoch nsat ',iepoch,nsat 
         call readc4 ( ifile
     .,                msat,mprn
     .,                itepch,iyr,idoy,sod,rclock
     .,                okmet,zendel,pres,temp,relhum,atmlod4
     .,                sitecd,latr_sph,lonr,radius
     .,                l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,                nsave,save  )
                                                           
         do i=1,3
           atmload(i,istat) = dble(atmlod4(i))
         enddo

         do i = 1, nsat
           ierrfl(i) = 1
         enddo
         if(debug) print *,'GETCDA 2 A ',iepoch,(a(ii),ii=1,3)
 

         do 30 jsat = 1, msat
            do 20 i = 1, nsat
               if (isprn(i).eq.mprn(jsat)) then
                  call readc5 ( ifile 
     .                        , iprn 
     .                        , el(i),elvdot 
     .                        , azim(i),azmdot,nadang
     .                        , atmdel 
     .                        , svcepc,svcl1 
     .                        , tau,drate 
     .                        , ierrfl(i),data_flag(i) 
     .                        , ndat 
     .                        , data(1,i) 
     .                        , omc(1,i) 
     .                        , isnr(1,i),ampl1(i),ampl2(i) 
     .                        , nspare,spare 
     .                        , npartc 
     .                        , tpart(1,istat,i) )     
                                          
        if(debug) print *,'GETCDA 3 A ',iepoch,(a(ii),ii=1,3)
 
c** zero out the L2 values for single-frequency receivers
        if( lwave(istat,i,2).eq.0 ) then
            data(2,i) = 0.d0
            data(4,i) = 0.d0
            omc(2,i) = 0.d0
            omc(4,i) = 0.d0 
         endif
                 
c** rwk 070328: 'iepoch' is passed to this routine only for debugging.  To
c               avoid a compiler warning for non-use, add this meaningless statement:
         j= iepoch       
      if(debug) print *,'GETCDA iepoch stat sat omc el ierrfl data_flag'
     .        ,iepoch,istat,i,omc(1,i),el(i),ierrfl(i),data_flag(i)

c***      
cd      print *,'iepoch svcepc svcl1 tau ierrfl omc '
cd     .      , iepoch,svcepc,svcl1,tau,ierrfl(i),omc(1,i)
cd      if( iepoch.eq.40 ) stop
cd        compute the atmospheric gradient partial
                  call atm_grad_part (el(i),azim(i),partial)
cd                  print*,'iepoch,istat,i,el(i),azim(i),partial: ',
cd     .                    iepoch,istat,i,el(i),azim(i),partial
                  do j = 1,2
                    gpart(istat,i,j) = partial(j)
                  enddo
                  goto 30
               endif
          if(debug) print *,'GETCDA gpart ',(gpart(istat,i,j),j=1,2)
 20         continue
 30      continue
         if(debug) print *,'GETCDA 4 A ',iepoch,(a(ii),ii=1,3)
 

      return
      end
