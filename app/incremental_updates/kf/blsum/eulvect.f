      program eulvect

c  program to calculate the Euler vectors for a series of plates given input information
c  of site coordinates (need only be approximate), site velocities and their uncertainties.
c  This program is designed to replace the DCCM program of Chuck DeMets and is based on
c  the code of TAH from the glorg/est_plate subroutine.
c
c  This program will calculate the absolute Euler vectors and also the relative vectors
c  between all plates. It will print some info to the screen and will also write out
c  to a specified file the Euler vectors and the VCV (in spherical coords) of the vectors.
c  The latter information is required by program EULER to calculate the uncertainties of
c  predicted velocities using the Euler vectors computed here.
c
c  The input file for this program has the following format:
c
c  Line 1:  Header line
c  Line 2:  Plates (a20 then free format) <--  the following lines are a priori values
c  Line 3:  "Plate Name"   "Plate Num"     "lat"  "long"    "omega"
c    .
c    .
c    .
c  Line X:  end
c           Site    Plate Num      X   Y    Z     Vx    Vy    Vz    Sx   Sy   Sz   Sxy   Sxz   Syz
c            .
c            .
c            .
c            end
c
c  P. Tregoning
c  13 October 2000
c
c  Mods:
c 
c  PT010912: increase maxplate from 10 to 20  and maxsit from 100 to 200
c  PT020304: decrease to 5 and 50 so it runs on amzer ...
c  SCM2010-10-29: modified to compile in blsum. Added some features...

      implicit none

      integer maxplate,maxsit
      parameter (maxplate = 50,maxsit = 1500)
      integer i,ierr,j,k,platenum(maxplate),nplates,nsit,psit
     .        ,sitnum(maxsit,2),uin,uout

      real*8 xyz(maxsit,3),velxyz(maxsit,3)
     .       ,covxyz(maxsit,3,3),ev_apr(maxplate,3)
     .       ,amat(maxsit*3,maxplate*3),bmat(maxsit*3)
     .       ,wmat(maxsit*3,maxsit*3)
     .       ,amatt(maxplate*3,maxsit*3),tmp1(maxsit*3,maxsit*3)
     .       ,tmp2(maxsit*3,maxsit*3),tmp3(maxplate*3)  
     .       ,tmpev_xyz(maxplate*3),ev_sph(maxplate,3)
     .       ,lat,long,w,pi
     .       ,ev_xyz(maxplate,3),covsoln(maxplate*3,maxplate*3)
     .       ,covsph(maxplate*3,maxplate*3)
c     .       ,relev_xyz(maxplate*maxplate,3)
c     .       ,relev_sph(maxplate*maxplate,3) 
c     .       ,cov
c
      character infile*30,outfile*20,ent*1,line*150,plates(maxplate)*20
     .          ,site(maxsit)*10

c   define unit numbers for the two files
      uin = 10
      uout = 20

c  the usual stuff first - get the filenames and open the files
      call getarg(1,infile)
      if(infile(1:1).eq.' ')then
        print*,' Runstring: eulvect <infile> <outfile>'
        stop 'Incorrect Runstring. Program stopped.'
      endif
    
      open(unit=uin,file=infile,status='old',iostat=ierr)
      if(ierr.ne.0)then
        print*,' Error opening file ',infile
        stop 'Program stopped'
      endif

      call getarg(2,outfile)
      open(unit=uout,file=outfile,status='unknown',iostat=ierr)
      if(ierr.ne.0)then
        print*,' Error opening file ',outfile
        stop 'Program stopped'
      endif
  
c  do a quick read of the input file to determine how many plates and sites there are.
c  I do this so that when I store the info in the next read I can actually dimension
c  the arrays with the number of sites and plates (and therefore can use printmat to
c  view the contents
      call read_input(uin,uout,maxsit,maxplate,plates,platenum,ev_apr
     .               ,xyz,velxyz,covxyz,site,nsit,psit,sitnum,nplates)
      rewind(uin)

c now read the information in the input file and store it in the arrays
      call read_input(uin,uout,nsit,nplates,plates,platenum,ev_apr
     .               ,xyz,velxyz,covxyz,site,nsit,psit,sitnum,nplates)

c  now form up the partials and the like
      call formup(uin,uout,plates,platenum,ev_apr
     .               ,xyz,velxyz,covxyz,site,nsit,psit,sitnum,nplates
     .               ,amat,bmat,wmat,tmp1)

*               T -1       T -1
****  Now form A W  A and A W  b
*

c  transpose amat
      call transp(amat,amatt,psit*3,nplates*3)

      call matmult(amatt,wmat,tmp1,nplates*3,psit*3,psit*3)

      call matmult(tmp1,amat,tmp2,nplates*3,psit*3,nplates*3)


c invert this (this gives the a posteriori covariance matrix of the parameters!)
      call invert(tmp2,covsoln,nplates*3,nplates*3)
c      call printmat(covsoln,nplates*3,nplates*3,'covsoln_in_main ')

c now the second term
      call matmult(tmp1,bmat,tmp3,nplates*3,psit*3,1)

c eh bah voila - la solution! 
      call matmult(covsoln,tmp3,tmpev_xyz,nplates*3,nplates*3,1)

c update the solution
      call solnupd(psit,nplates,ev_apr,tmpev_xyz,ev_xyz  
     .              ,plates,ev_sph,covsoln,covsph)

c and outout the soln. The relative poles and VCV are also calculated here and
c written out one at a time 
      call writeout(uout,psit,nplates,ev_xyz
     .              ,plates,ev_sph,covsoln,covsph)

c  now calculate the a posteriori residuals
      call apost(nsit,psit,nplates,plates,platenum,ev_sph,xyz,velxyz,
     .           site,sitnum,covxyz,uout)

c  calculate site velocities in all plate frames
      call pvel(nsit,nplates,plates,platenum,ev_sph,xyz,velxyz,site
     .           ,sitnum,covxyz)

c'est tout!
      end

      subroutine read_input(uin,uout,maxsit,maxplate,plates,platenum
     .                      ,ev_apr,xyz,velxyz,covvelxyz,site
     .                      ,nsit,psit,sitnum,nplates)

c  subroutine to read the input data and store it away in arrays
c
c  IN:  
c      uin             : unit number of input file
c      uout            : unit number of output file
c      maxsit,maxplate : maximum numbers of each
c
c  OUT:
c      plates          : Plate names 
c      platenum        : plate number
c      ev_apr          : a priori euler vectors for each plate (in spherical coords) (deg and deg/yr)
c      xyz             : site coordinates    (in metres)
c      velxyz          : site velocities     (m/yr)
c      covvelxyz       : velocity covariance matrices  (m**2/yr)
c
c  P. Tregoning 
c  13 October 2000
c
c  MODS:
c  001013: this routine is called twice, the first time to find out how many
c          plates and sites there are (when maxsit and maxplate are the actual
c          dimensioned values in the main program). 
c
c          The second time these values are actually nsit and nplate; therefore
c          the data stored in the arrays can be viewed using printmat in any
c          subroutine where the arrays are passed through and dimensioned using
c          nsit and nplates

      integer uin,uout,maxplate,maxsit,ierr,i,j,k,platenum(maxplate)
     .        ,nplates,nsit,psit,pltmp,sitnum(maxsit,2)
     .        ,count

      real*8 xyz(maxsit,3),velxyz(maxsit,3)
     .       ,covvelxyz(maxsit,3,3),ev_apr(maxplate,3),covtmp(6)
      character infile*20,outfile*20,ent*1,line*180,plates(maxplate)*20
     .          ,site(maxsit)*10

c  The input file for this program has the following format:
c
c  Line 1:  Header line
c  Line 2:  Plates (a20 then free format) <--  the following lines are a priori values
c  Line 3:  "Plate Name"   "Plate Num"     "lat"  "long"    "omega"
c    .
c    .
c    .
c  Line X:  end
c           Site    Plate Num      X   Y    Z     Vx    Vy    Vz    Sx   Sy   Sz   Sxy   Sxz   Syz
c            .
c            .
c            .
c            end

      nplates = 0
      nsit = 0
      psit = 0 
      line = ' '
      count = 0

c  OK, read the input header line and write it straight to the output file
      read(uin,10)line
10    format(a)
      write(uout,10)line

c  next line can be ignored
      read(uin,10)line
      
c now loop through reading the a priori Euler vectors for the plates until "end" is read  
      do while(line(1:3).ne.'end')
        read(uin,10)line
        if(line(1:3).ne.'end')then  
          nplates = nplates + 1                     
          read(line(1:20),'(a)')plates(nplates)  
          read(line(21:80),*)platenum(nplates)
     .                ,ev_apr(nplates,1)
     .                ,ev_apr(nplates,2),ev_apr(nplates,3)  

c change deg/My to deg/yr  
c  No, don't!
c yes do!!
          ev_apr(nplates,3) = ev_apr(nplates,3) *1.d-6   

c        print*,plates(nplates),platenum(nplates)
c     .                ,ev_apr(nplates,1)
c     .                ,ev_apr(nplates,2),ev_apr(nplates,3)

        endif
      enddo

c now loop through the list of site entries and store the coordinate/velocity information    
      line = ' '
      do while(line(1:3).ne.'end')
        read(uin,10)line
        if(line(1:3).ne.'end')then  
          read(line(11:20),*)pltmp        

c check that the plate number exists in the list
          if(pltmp.gt.nplates)then
c            print*,'Site not assigned to a plate: ', line 
            nsit = nsit + 1  
            sitnum(nsit,1) = nsit
            sitnum(nsit,2) = pltmp
	  else
            nsit = nsit + 1  
            sitnum(nsit,1) = nsit
            sitnum(nsit,2) = pltmp
	    psit = psit + 1
          endif

c  now that we know which plate it is on, store the values in the other arrays
          read(line(1:10),'(a)')site(nsit)
          read(line(11:180),*)pltmp,(xyz(nsit,j),j=1,3)
     .                ,(velxyz(nsit,j),j=1,3),(covtmp(j),j=1,6)        

c  now store the covariance array in its proper form  
c  The values in the input file are listed as
c   Sx  Sy  Sz  Sxy  Sxz Syz
          do i=1,3
            covvelxyz(nsit,i,i) = covtmp(i)**2
          enddo


c  PT001019: The covariance off-diagonal terms listed in the input file are actually
c            "sigmas" rather than "variance" terms as calculated from 
c
c   omega(XY) = omega**2(XY)/(dsqrt(var(X)*var(Y))
c
c            I need to undo this. 
c
c PT010816: The actual off-diagonal covariance terms are now output to the glorg file
c           so I no longer need to do this stuff! Just copy them from covtmp straight
c           into covvelxyz

c          covvelxyz(nsit,1,2) = covtmp(4)
c     .                *dsqrt(covvelxyz(nsit,1,1)*covvelxyz(nsit,2,2))
c          covvelxyz(nsit,1,3) = covtmp(5)
c     .                *dsqrt(covvelxyz(nsit,1,1)*covvelxyz(nsit,3,3))
c          covvelxyz(nsit,2,3) = covtmp(6)
c     .                *dsqrt(covvelxyz(nsit,2,2)*covvelxyz(nsit,3,3))
                                    
          covvelxyz(nsit,1,2) = covtmp(4)
          covvelxyz(nsit,1,3) = covtmp(5)
          covvelxyz(nsit,2,3) = covtmp(6)

c  and fill the upper-diagonal space
          covvelxyz(nsit,2,1) = covvelxyz(nsit,1,2)
          covvelxyz(nsit,3,1) = covvelxyz(nsit,1,3)
          covvelxyz(nsit,3,2) = covvelxyz(nsit,2,3) 

c          print*,' covariance matrix for ',site(nsit)
c          print*,(covvelxyz(nsit,1,j),j=1,3)
c          print*,(covvelxyz(nsit,2,j),j=1,3)
c          print*,(covvelxyz(nsit,3,j),j=1,3)
        endif
      enddo

      return
      end

      subroutine writeout(uout,psit,nplates,ev_xyz
     .                    ,plates,ev_sph,covsoln,covsph)

c  subroutine to output the absolute Euler vectors and their VCV matrices in
c  spherical coordinates. The routine also calculates the relative poles for
c  each possible plate combination and then outputs this in spherical coords.
c
c  P. Tregoning
c  13 October 2000

      implicit none

      integer uout,psit,nplates,ierr,i,j,k,pl1,pl2
      real*8 ev_xyz(nplates,3)
     .      ,ev_sph(nplates,3)
     .      ,covsoln(nplates*3,nplates*3),covsph(nplates*3,nplates*3) 
     .      ,lat,long,w,w2,pi,rt(3,3)
     .      ,rotmat1(3,3),relxyz(3),relsph(3),tmp1(6,6),tmpcov(3,3)
     .      ,rotmat2(3,3),smaj,smin,az,tmpcov2(6,6),jac(3,6),jact(6,3)
     .      ,relcov(3,3)    
      character plates(nplates)*20

      pi=4.d0*datan(1.d0)


c  form the jacobian matrix for calculating the VCV of relative poles
      do i=1,3
        do j=1,3
          jac(i,j) = 0.d0
          jac(i,j+3) = 0.d0 
          if(i.eq.j)then
            jac(i,j) = 1.d0
            jac(i,j+3) = -1.d0
          endif
        enddo
      enddo
      
      write(*,70)
      write(uout,70)
70    format('***********************************************',/,
     .       ' A B S O L U T E    V E C T O R S',/,
     .       '***********************************************')
      
c  now loop through and output all the Absolute Poles
      do i=1,nplates     

c extract the spherical VCV from the  matrix of spherical VCV's
        do j=1,3  
          do k=1,3
            tmpcov(j,k) =  covsph(i*3-3+j,k)
          enddo
        enddo

c  convert it to semi-maj and semi-min error ellipse
        call errell(tmpcov(1,1),tmpcov(2,2),tmpcov(1,2),smaj,smin,az)

c  now write out the plate name and Euler vector along with the smaj, smin and azimuth
c and output the relative Euler vectors
        write(*,95)'Plate','Lat','Long','Rate',
     .              'Semi Major','Semi Minor','Azim','Rate Unc'
        write(uout,95)'Plate','Lat','Long','Rate',
     .                'Semi Major','Semi Minor','Azim','Rate Unc'
95      format(1x,a5,15x,5(a8,1x),a7,1x,a9)
        write(*,96)'deg','deg',' deg/Myr',
     .                'deg','deg','deg','deg/Myr'
        write(uout,96)'deg','deg',' deg/Myr',
     .                'deg','deg','deg','deg/Myr'
96      format(21x,5(a8,1x),a7,1x,a8)        
        write(*,100)plates(i),(ev_sph(i,j),j=1,3)
     .          ,smaj*180.d0/pi,smin*180.d0/pi,az
     .          ,dsqrt(tmpcov(3,3))*1.d6*180.d0/pi
        write(uout,100)plates(i),(ev_sph(i,j),j=1,3)
     .          ,smaj*180.d0/pi,smin*180.d0/pi,az
     .          ,dsqrt(tmpcov(3,3))*1.d6*180.d0/pi   
100     format(1x,a20,3(f8.3,1x),2(f8.2,1x),f7.1,1x,f8.4)

c and now the VCV in spherical coordinates   
c PT001019: this is the standard deviations, I need the variances
c        write(*,110)
cc row1
c     .  dsqrt(tmpcov(1,1))*180.d0/pi
c     .  ,tmpcov(1,2)/dsqrt(tmpcov(1,1)*tmpcov(2,2))
c     .  ,tmpcov(1,3)/dsqrt(tmpcov(1,1)*tmpcov(3,3)*1.d6)
cc row2
c     .  ,tmpcov(2,1)/dsqrt(tmpcov(1,1)*tmpcov(2,2))  
c     .  ,dsqrt(tmpcov(2,2))*180.d0/pi
c     .  ,tmpcov(2,3)/dsqrt(tmpcov(2,2)*tmpcov(3,3)*1.d6)
cc row3
c     .  ,tmpcov(3,1)/dsqrt(tmpcov(1,1)*tmpcov(3,3)*1.d6)
c     .  ,tmpcov(3,2)/dsqrt(tmpcov(2,2)*tmpcov(3,3)*1.d6)
c     .  ,dsqrt(tmpcov(3,3))*180.d6/pi

        write(*,125)'VCV in spherical coordinates'
        write(*,126)'Lat ','Long','Rate'
        write(*,110)
c row1
     .  tmpcov(1,1)
     .  ,tmpcov(1,2)
     .  ,tmpcov(1,3)
c row2
     .  ,tmpcov(2,1)
     .  ,tmpcov(2,2)
     .  ,tmpcov(2,3)
c row3
     .  ,tmpcov(3,1)
     .  ,tmpcov(3,2)
     .  ,tmpcov(3,3)
     	
        write(uout,125)'VCV in spherical coordinates'
        write(uout,126)'Lat ','Long','Rate'
        write(uout,110)tmpcov(1,1),tmpcov(1,2),tmpcov(1,3),tmpcov(2,1),
     .       tmpcov(2,2),tmpcov(2,3),tmpcov(3,1),tmpcov(3,2),tmpcov(3,3)
     
110     format(3(3e19.9,/))

      enddo

      write(*,112)
      write(uout,112)
112   format('***********************************************',/,
     .       ' R E L A T I V E    V E C T O R S',/,
     .       '***********************************************')

c  now calculate all the relative poles
      do pl1=1,nplates-1
        do pl2=pl1+1,nplates
c calculate the motion of the second plate wrt the first
c subtract the xyz of the second plate from the first plate
          do k=1,3
            relxyz(k) = ev_xyz(pl1,k) - ev_xyz(pl2,k)
          enddo

c  now convert to spherical
          call carsph(relxyz,relsph,1,rotmat1,rotmat2)
          relsph(1) = relsph(1)*180.d0/pi
          relsph(2) = relsph(2)*180.d0/pi
          relsph(3) = relsph(3)*180.d6/pi

c  now compute the VCV of the spherical representation. First extract the 3x3 vector
c  for each pole and store them in a 6 x 6 array
          do j=1,3  
            do k=1,3
              tmpcov2(j,k) = covsoln(pl1*3 - 3 + j,pl1*3 - 3 + k) 
              tmpcov2(j+3,k+3) = covsoln(pl2*3 - 3 + j,pl2*3 - 3 + k) 
            enddo
          enddo  
c                                                 
c                                                t
c  multiply it out   VCV     =  jac   VCV     Jac    where Jac = jacobian for subtracting two poles
c                       rel               xyz     
  
          call transp(jac,jact,3,6)
          call matmult(jac,tmpcov2,tmp1,3,6,6)
          call matmult(tmp1,jact,relcov,3,6,3)
                                                      

c  finally, convert this to a VCV of spherical coordinates
c                                                 -1
c                                  -1            t
c  multiply it out   VCV     =  jac   VCV     Jac    where Jac = jacobian SPH -> XYZ
c                       sph               xyz     
  
          call transp(rotmat2,rt,3,3)
          call matmult(relcov,rt,tmp1,3,3,3)
          call matmult(rotmat2,tmp1,relcov,3,3,3)

c  convert it to semi-maj and semi-min error ellipse
          call errell(relcov(1,1),relcov(2,2),relcov(1,2),smaj,smin,az)


c and output the relative Euler vectors
          write(*,115)'Plate 1','Plate 2','Lat','Long','Rate',
     .                'Semi Major','Semi Minor','Azim','Rate Unc'
          write(uout,115)'Plate 1','Plate 2','Lat','Long','Rate',
     .                'Semi Major','Semi Minor','Azim','Rate Unc'
115       format(1x,a7,14x,a7,14x,5(a8,1x),a7,1x,a8)
          write(*,116)'deg','deg',' deg/Myr',
     .                'deg','deg','deg','deg/Myr'
          write(uout,116)'deg','deg',' deg/Myr',
     .                'deg','deg','deg','deg/Myr'
116       format(43x,5(a8,1x),a7,1x,a8)
          write(*,120)plates(pl2),plates(pl1),relsph
     .          ,smaj*180.d0/pi,smin*180.d0/pi,az
     .          ,dsqrt(relcov(3,3))*1.d6*180.d0/pi
          write(uout,120)plates(pl2),plates(pl1),relsph
     .          ,smaj*180.d0/pi,smin*180.d0/pi,az
     .          ,dsqrt(relcov(3,3))*1.d6*180.d0/pi
120       format(1x,a20,1x,a20,1x,3(f8.3,1x),2(f8.2,1x),f7.1,1x,f8.4)

c and now the VCV in spherical coordinates   
          write(*,125)'VCV in spherical coordinates'
          write(uout,125)'VCV in spherical coordinates'
125       format(/,1x,a)
          write(*,126)'Lat ','Long','Rate'
          write(uout,126)'Lat ','Long','Rate'
126       format(1x,3(a10,9x))
	  
          write(*,130)
cc row1
c     .    dsqrt(relcov(1,1))*180.d0/pi
c     .    ,relcov(1,2)/dsqrt(relcov(1,1)*relcov(2,2))
c     .    ,relcov(1,3)/dsqrt(relcov(1,1)*relcov(3,3)*1.d6)
cc row2
c     .    ,relcov(2,1)/dsqrt(relcov(1,1)*relcov(2,2))  
c     .    ,dsqrt(relcov(2,2))*180.d0/pi
c     .    ,relcov(2,3)/dsqrt(relcov(2,2)*relcov(3,3)*1.d6)
cc row3
c     .    ,relcov(3,1)/dsqrt(relcov(1,1)*relcov(3,3)*1.d6)
c     .    ,relcov(3,2)/dsqrt(relcov(2,2)*relcov(3,3)*1.d6)
c     .    ,dsqrt(relcov(3,3))*180.d6/pi
c130       format(3(3f14.9,/))
c row1
     .    relcov(1,1)
     .    ,relcov(1,2)
     .    ,relcov(1,3)
c row2
     .    ,relcov(2,1)
     .    ,relcov(2,2)
     .    ,relcov(2,3)
c row3
     .    ,relcov(3,1)
     .    ,relcov(3,2)
     .    ,relcov(3,3)
     
          write(uout,130)relcov(1,1),relcov(1,2),relcov(1,3),relcov(2,1)
     .      ,relcov(2,2),relcov(2,3),relcov(3,1),relcov(3,2),relcov(3,3)

130       format(3(3e19.9,/))
        enddo
      enddo

      return
      end

        subroutine transp(R,Rt,m,n)

c  transposes a mxn matrix
c  Written by Paul Tregoning
c  10/6/92

        implicit none
        integer i,j,m,n
        real*8 R(m,n),Rt(n,m)

c  set Rt to zero initially
        do 5 i=1,n
          do 8 j=1,m
             Rt(i,j)=0.0
8         continue
5       continue

c  perform the transpose
        do 10 i=1,m
        do 20 j=1,n
            Rt(j,i) = R(i,j)
20        continue
10      continue

        return
        end

        subroutine matmult(M1,M2,M3,l,m,n)

c  multiplies lxm matrix M1 by a mxn matrix M2 to give a lxn matrix M3
c  written by Paul Tregoning
c  9th June 1992

        implicit none
        integer i,j,k,l,m,n
        real*8 M1(l,m),M2(m,n),M3(l,n),temp

        temp = 0.d0

        do 10 i=1,l
          do 30 k=1,n
            do 20 j=1,m
              temp = M2(j,k) * M1(i,j) + temp  
c              print*,i,k,j,M2(j,k) * M1(i,j),temp
20          continue

          M3(i,k) = temp
          temp = 0.d0
30       continue
10      continue


        return
        end

      subroutine apost(nsit,psit,nplates,plates,platenum,ev_sph,xyz,
     .                 velxyz,site,sitnum,covxyz,uout)

c  subroutine to calculate the a posteriori residuals and all the
c  statistics associated with LS solutions. Variance factors, chi**2/dof etc
c
c  P. Tregoning
c  17 october 2000

      implicit none

      integer i,j,isit,iplate,nsit,psit,nplates,platenum(nplates)
     .       ,sitnum(nsit,2),pnum,nperplate,uout
      real*8 ev_sph(nplates,3),xyz(nsit,3),velxyz(nsit,3)
     .       ,cvelxyz(3),tmpxyz(3),loc_coord(3),rotmat(3,3)
     .       ,covxyz(nsit,3,3),tmpcov(3,3),tmpcovneu(3,3)
     .       ,tmpvelxyz(3),velneu(3),residxyz(3),tmp1(3,3),rt(3,3)
     .       ,pi,rdeg,cvelneu(3),ovelneu(3),resid(2),chi2_tot
     .       ,wrms_n,wrms_e,sumw_n,sumw_e,chi2_n,chi2_e
     .       ,wrms_tot

      character site(nsit)*10,plates(nplates)*20

      pi=4.d0*datan(1.d0)
      rdeg = pi/180.d0

      write(*,85)
      write(uout,85)
85    format('Plate Fit Statistics and Residuals',/)

c  set up a loop through all the plates. Calculate the model velocities for
c  each site on that plate and subtract it from the input velocity
      do iplate = 1,nplates 
      
c initialize stats variables      
        nperplate = 0
        resid(1) = 0.d0
        resid(2) = 0.d0
        wrms_e = 0.d0
        wrms_n = 0.d0
        sumw_n = 0.d0
        sumw_e = 0.d0
	
c  write a header line for the residuals output
        write(*,90)'Site','Plate','dVn','dVe','dVu','Mod Vn','Mod Ve'
     .      ,'Mod Vu','Obs Vn','Obs Ve','Obs Vu','sVn','sVe','sVu'
        write(uout,90)'Site','Plate','dVn','dVe','dVu','Mod Vn','Mod Ve'
     .      ,'Mod Vu','Obs Vn','Obs Ve','Obs Vu','sVn','sVe','sVu'
90      format(3x,a4,10x,a5,10x,3(3x,a3,2x),3(2x,a6),3(2x,a6)
     .      ,2x,3(a3,3x)) 
     
        do isit = 1,nsit  
          if(sitnum(isit,2).eq.platenum(iplate))then 
            nperplate = nperplate + 1
            pnum = platenum(iplate)
c it is on this plate
            call calcvel(xyz,nplates,pnum,nsit,isit,ev_sph(pnum,1)
     .        ,ev_sph(pnum,2),ev_sph(pnum,3)*1.d-6
     .        ,cvelxyz)

c  now calculate the residual (and store the xyz coords)
            do i = 1,3  
c              print*,' resid ',i,velxyz(isit,i),cvelxyz(i)
              residxyz(i) = velxyz(isit,i) - cvelxyz(i) 
              tmpxyz(i) = xyz(isit,i) 
              tmpvelxyz(i) = velxyz(isit,i)
            enddo
c  and convert it to NEU  
c PT010816: try using XYZ_to_GEOD instead      
c            call XYZ_to_NEU(rotmat,tmpxyz,loc_coord)
c        call printmat(rotmat,3,3,'rotmat_NEU ')
            call XYZ_to_GEOD(rotmat,tmpxyz,loc_coord)
            call matmult(rotmat,residxyz,velneu,3,3,1)

c  also convert the VCV     to VCV     for the input velocity
c                      xyz        neu    
c
            do i = 1,3
              do j=1,3
                tmpcov(i,j) = covxyz(isit,i,j) 
                tmp1(i,j) = 0.d0
              enddo
            enddo  

            call transp(rotmat,rt,3,3)
            call matmult(tmpcov,rt,tmp1,3,3,3)
            call matmult(rotmat,tmp1,tmpcovneu,3,3,3)

c  convert the calculate velocity to NEU to check it's right!
            call matmult(rotmat,cvelxyz,cvelneu,3,3,1)

c  convert the original velocity to NEU
            call matmult(rotmat,tmpvelxyz,ovelneu,3,3,1)

            write(*,100)site(isit),plates(platenum(iplate))
     .        ,(velneu(j)*1.d3,j=1,3)
     .        ,(cvelneu(j)*1.d3,j=1,3)
     .        ,(ovelneu(j)*1.d3,j=1,3)
     .        ,((dsqrt(tmpcovneu(j,j))*1.d3),j=1,3)
            write(uout,100)site(isit),plates(platenum(iplate))
     .        ,(velneu(j)*1.d3,j=1,3)
     .        ,(cvelneu(j)*1.d3,j=1,3)
     .        ,(ovelneu(j)*1.d3,j=1,3)
     .        ,((dsqrt(tmpcovneu(j,j))*1.d3),j=1,3)
100         format(a,1x,a,3f8.2,6f8.2,3f6.1)
c            stop ' stopped in apost'

c  now sum the residuals for north and east. Ignore the U component completely.  
c  PT011026: also compute the rms for N and E. 
c  That is WRMS = Sum( data[i] * data[i] * weights[i]) / Sum( data[i] * weights[i] )
c
            resid(1) = resid(1) + (velneu(1)**2/tmpcovneu(1,1))
            resid(2) = resid(2) + (velneu(2)**2/tmpcovneu(2,2)) 
            sumw_n = sumw_n + 1/tmpcovneu(1,1)
            sumw_e = sumw_e + 1/tmpcovneu(2,2)
          endif   
        enddo

c  PT011029: wrms = sqrt [sum ( V(i)*V(i)/W(i)) /n ]
          chi2_tot = (resid(1) + resid(2))/(nperplate*2.d0-3.d0)

          wrms_n = dsqrt(nperplate/sumw_n*(resid(1)/(nperplate)))*1.e3
          wrms_e = dsqrt(nperplate/sumw_e*(resid(2)/(nperplate)))*1.e3
          wrms_tot = dsqrt(nperplate/sumw_n*(resid(1)/(nperplate)) +
     .                nperplate/sumw_e*(resid(2)/(nperplate)))*1.e3

          write(*,110)'Chi**2 for ',plates(iplate),chi2_tot
     .      ,' : DoF = ',nperplate*2-3, ' : Chi**2/DoF = '
     .      ,Chi2_tot/(nperplate*2.d0-3.d0)
          write(uout,110)'Chi**2 for ',plates(iplate),chi2_tot
     .      ,' : DoF = ',nperplate*2-3, ' : Chi**2/DoF = '
     .      ,Chi2_tot/(nperplate*2.d0-3.d0)

110       format(a,a,1x,f7.2,a,i3,a,f7.2)

          write(*,111)'wrms N, E  ',plates(iplate),'   '
     .       ,wrms_n,wrms_e,' mm/yr. Total wrms is ',wrms_tot,' mm/yr'
          write(uout,111)'wrms N, E  ',plates(iplate),'   '
     .       ,wrms_n,wrms_e,' mm/yr. Total wrms is ',wrms_tot,' mm/yr'
111       format(a,a,a,2(f5.2,1x),a,f5.2,a,/)

      enddo

      return
      end
      
      subroutine pvel(nsit,nplates,plates,platenum,ev_sph,xyz,velxyz
     .                 ,site,sitnum,covxyz)

c  subroutine to calculate site velocities in all estimated plate frames
c
c S. McClusky
c  17 october 2010

      implicit none

      integer i,j,isit,iplate,nsit,nplates,platenum(nplates)
     .       ,sitnum(nsit,2),pnum,uout,trimlen,ierr
      real*8 ev_sph(nplates,3),xyz(nsit,3),velxyz(nsit,3)
     .       ,cvelxyz(3),tmpxyz(3),loc_coord(3),rotmat(3,3)
     .       ,covxyz(nsit,3,3),tmpcov(3,3),tmpcovneu(3,3)
     .       ,tmpvelxyz(3),velneu(3),residxyz(3),tmp1(3,3),rt(3,3)
     .       ,cvelneu(3),ovelneu(3),pi,rdeg

      character site(nsit)*10,plates(nplates)*20,outfile*20,pname*20

      pi=4.d0*datan(1.d0)
      rdeg = pi/180.d0
      uout = 30
     
c  write a header line for the residuals output
c      write(*,90)'Site','Plate','dVn','dVe','dVu','Mod Vn','Mod Ve'
c     .      ,'Mod Vu','Obs Vn','Obs Ve','Obs Vu','sVn','sVe','sVu'
c90    format(/,3x,a4,10x,a5,10x,3(3x,a3,2x),3(2x,a6),3(2x,a6)
c     .      ,2x,3(a3,3x))    

c  set up a loop through all the plates. Calculate the model velocities for
c  each site on that plate and subtract it from the input velocity
      print*,' '  
      do iplate = 1,nplates 
	pname = plates(platenum(iplate))
        outfile = pname(1:trimlen(pname))//".vel"
        print*,'Writing: ',outfile      
        open(unit=uout,file=outfile,status='unknown',iostat=ierr)
        if(ierr.ne.0)then
          print*,' Error opening file ',outfile
          stop 'Program stopped'
        endif
        write(uout,250)plates(platenum(iplate))
 250    format('*  ',a8)
        write(uout,260)
 260    format('*  Long         Lat        Evel    Nvel    dEv  ',
     .         '   dNv    E +-    N +-    Rne      Hvel     dHv    ',
     .         'H +-  Site',/,
     .         '*  deg          deg       mm/yr   mm/yr   mm/yr   ',
     .         'mm/yr   mm/yr   mm/yr            mm/yr   mm/yr',
     .         '   mm/yr')
        do isit = 1,nsit  
            pnum = platenum(iplate)

            call calcvel(xyz,nplates,pnum,nsit,isit,ev_sph(pnum,1)
     .        ,ev_sph(pnum,2),ev_sph(pnum,3)*1.d-6
     .        ,cvelxyz)

c  now calculate the residual (and store the xyz coords)
            do i = 1,3  
c              print*,' resid ',i,velxyz(isit,i),cvelxyz(i)
              residxyz(i) = velxyz(isit,i) - cvelxyz(i) 
              tmpxyz(i) = xyz(isit,i) 
              tmpvelxyz(i) = velxyz(isit,i)
            enddo
	    
c  and convert it to NEU  
c        call printmat(rotmat,3,3,'rotmat_NEU ')
            call XYZ_to_GEOD(rotmat,tmpxyz,loc_coord)
            call matmult(rotmat,residxyz,velneu,3,3,1)

c  also convert the VCV     to VCV     for the input velocity
c                      xyz        neu    
c
            do i = 1,3
              do j=1,3
                tmpcov(i,j) = covxyz(isit,i,j) 
                tmp1(i,j) = 0.d0
              enddo
            enddo  

            call transp(rotmat,rt,3,3)
            call matmult(tmpcov,rt,tmp1,3,3,3)
            call matmult(rotmat,tmp1,tmpcovneu,3,3,3)

c  convert the calculate velocity to NEU to check it's right!
            call matmult(rotmat,cvelxyz,cvelneu,3,3,1)

c  convert the original velocity to NEU
            call matmult(rotmat,tmpvelxyz,ovelneu,3,3,1)

c  output site velocities rotated into each plate frame
          if(sitnum(isit,2).eq.platenum(iplate))then 
              write(uout,500) loc_coord(2)*180/pi, 
     .        (pi/2-loc_coord(1))*180/pi,
     .        velneu(2)*1.d3,velneu(1)*1.d3,
     .        ovelneu(2)*1.d3,ovelneu(1)*1.d3,
     .        dsqrt(tmpcovneu(2,2))*1.d3,dsqrt(tmpcovneu(1,1))*1.d3,
     .        tmpcovneu(1,2)/(dsqrt(tmpcovneu(1,1))*
     .        dsqrt(tmpcovneu(2,2))), 
     .        ovelneu(3)*1.d3,velneu(3)*1.d3,dsqrt(tmpcovneu(3,3))*1.d3,
     .        site(isit)
c              write(*,100)site(isit),plates(platenum(iplate))
c     .        ,(velneu(j)*1.d3,j=1,3)
c     .        ,(cvelneu(j)*1.d3,j=1,3)
c     .        ,(ovelneu(j)*1.d3,j=1,3)
c     .        ,((dsqrt(tmpcovneu(j,j))*1.d3),j=1,3)
          else
              write(uout,510) loc_coord(2)*180/pi, 
     .        (pi/2-loc_coord(1))*180/pi,
     .        velneu(2)*1.d3,velneu(1)*1.d3,
     .        ovelneu(2)*1.d3,ovelneu(1)*1.d3,
     .        dsqrt(tmpcovneu(2,2))*1.d3,dsqrt(tmpcovneu(1,1))*1.d3,
     .        tmpcovneu(1,2)/(dsqrt(tmpcovneu(1,1))*
     .        dsqrt(tmpcovneu(2,2))), 
     .        ovelneu(3)*1.d3,velneu(3)*1.d3,dsqrt(tmpcovneu(3,3))*1.d3,
     .        site(isit)     
c              write(*,110)site(isit),plates(platenum(iplate))
c     .        ,(velneu(j)*1.d3,j=1,3)
c     .        ,(cvelneu(j)*1.d3,j=1,3)
c     .        ,(ovelneu(j)*1.d3,j=1,3)
c     .        ,((dsqrt(tmpcovneu(j,j))*1.d3),j=1,3)
     
c  write globk velocity file format output
 500          format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8,'*')
 510          format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8,' ')
     
c  write eulvect format output
c 100         format(a,1x,a,3f8.2,6f8.2,3f6.1,' *')
c 110         format(a,1x,a,3f8.2,6f8.2,3f6.1)
          endif
     
        enddo
	close(uout)
      enddo

      return
      end

      subroutine formup(uin,uout,plates,platenum,ev_apr
     .              ,xyz,velxyz,covvelxyz,site,nsit,psit,sitnum,nplates
     .              ,amat,bmat,wmat,tmp)


c  subroutine to calculate the partials and form the matrices ready for
c  the least squares solution
c
c  IN:  
c      uin             : unit number of input file
c      uout            : unit number of output file
c      maxsit,maxplate : maximum numbers of each
c      plates          : Plate names 
c      platenum        : plate number
c      ev_apr          : a priori euler vectors for each plate (in spherical coords)
c      xyz             : site coordinates
c      velxyz          : site velocities
c      covvelxyz       : velocity covariance matrices
c
c  OUT: 
c      amat   : a matrix
c      wmat   : weight matrix
c      bmat   : o-c matrix
c
c  P. Tregoning 
c  13 October 2000
      implicit none

      integer uin,uout,nplates,nsit,ierr,i,j,k,l,platenum(nplates)
     .        ,pltmp,sitnum(nsit,2)
     .        ,count,ps,pl,psit

      real*8 xyz(nsit,3),velxyz(nsit,3)
     .       ,covvelxyz(nsit,3,3),ev_apr(nplates,3),covtmp(6)
     .       ,amat(psit*3,nplates*3),bmat(psit*3),wmat(psit*3,psit*3)
     .       ,cvelxyz(3),tmp(psit*3,psit*3)
      character infile*20,outfile*20,ent*1,line*150,plates(nplates)*20
     .          ,site(nsit)*10


c  first, view a few matrices to see whether the data are stored properly
c      call printmat(ev_apr,nplates,3,'ev_apr ')
c      call printmat(velxyz,nsit,3,'velxyz ')
c      call printmat(xyz,nsit,3,'xyz ')
c      call printmat(covvelxyz(1,1,1),3,3,'covvelxyz ')

c  now comes the guts of the whole program ..... forming up the A matrix partials
c  The parameters in this whole process are 3 x rotations of the cartesian axes
c  for each plate. Therefore, we will have (3 x nplates) columns in the A matrix
c  and (3 x nsit) rows
      l = 0 
      do i=1,nsit
        if(sitnum(i,2) .ne. 99) then 
	  l = l + 1
	  ps = 3*(l-1)	  
          do j=1,nplates
            pl = 3*(j-1)
            if(sitnum(i,2).eq.j) then 
	      amat(ps+1,pl+2) =  xyz(i,3)
              amat(ps+1,pl+3) = -xyz(i,2)
      
              amat(ps+2,pl+1) = -xyz(i,3)
              amat(ps+2,pl+3) =  xyz(i,1)

              amat(ps+3,pl+1) =  xyz(i,2)
              amat(ps+3,pl+2) = -xyz(i,1)
            endif
          enddo
	  
c  now move the covariance elements into the right place in the VCV
          do j=1,3 
            do k=1,3  
              wmat(ps+j,ps+k) = covvelxyz(i,j,k)
            enddo
          enddo        

c  the B matrix (O-C) is straightforward. For each site we can calculate the
c  a priori velocity using the a priori Euler vector, and then subtract this
c  from the velocities inputted through the input file (ie stored in xyzvel)

c  first, calculate the a priori velocity from the model
          call calcvel(xyz,nplates,sitnum(i,2),nsit,i
     .        ,ev_apr(sitnum(i,2),1)
     .        ,ev_apr(sitnum(i,2),2),ev_apr(sitnum(i,2),3),cvelxyz)

c  now subtract it from the "observed" velocity
          do j=1,3
            bmat(ps + j) = velxyz(i,j) - cvelxyz(j) 
c            print*,'O-C ',bmat(ps + j)
          enddo
	  
        endif
	
      enddo
c      call printmat(amat,psit*3,nplates*3,'amat ')

c  now invert the weight matrix
      call invert(wmat,tmp,psit*3,psit*3)
c  and reallocate it
      do i = 1,psit*3
        do j=1,psit*3
          wmat(i,j) = tmp(i,j)
        enddo
      enddo

      return
      end

      subroutine invert(a,b,n,m)

c  subroutine stolen from someone else

      implicit none

      integer n,m,i,j,k
      real*8 a(m,m),b(m,m),z

      do  i=1,n
        do  j=1,n
          b(i,j) = 0.0
        enddo
        b(i,i) = 1.0
      enddo

      do k = 1,n
        do i = 1,n
	       if(i.ne.k)then
	         z = a(i,k)/a(k,k)
            do j=1,n
              a(i,j) = a(i,j) - a(k,j)*z
              b(i,j) = b(I,j) - b(k,j)*z
            enddo
	       endif
	     enddo
	     z = a(k,k)
	     do j=1,n
	       a(k,j) = a(k,j)/z
	       b(k,j) = b(k,j)/z
	     enddo
	   enddo

      return
      end

      subroutine carsph(ev_xyz,ev_sph,dirn,rotmat1,rotmat2)

c  subroutine to convert between spherical and cartesian representations of an Euler vector.
c  It is either lat,lon,w  in spherical coords or Rx, Ry, Rz rotations about the cartesian
c  axes. Routine will also output the jacobian matrix to transform from VCV spherical
c  to VCV cartesian. The inverse of this can be used to go from VCV car to VCV sph.
c                  
c  ARGUMENTS:
c      ev_xyz  :   Euler vector rotations as cartesian rotations
c      ev_sph  :   spherical representation of the Euler vector
c      rotmat  :   jacobian matrix (sph to xyz)
c      dirn    :   1 (XYZ to SPH) or 2 (SPH to XYZ)

c  P. Tregoning
c  16 october 2000    
c
c  MODS:
c  PT010816: fix bug that invert.f was clobbering rotmat1. This matrix was never used
c            again but at least the bug is now fixed!

      implicit none

      integer i,j,dirn
      real*8 ev_xyz(3),ev_sph(3),rotmat1(3,3),rotmat2(3,3),lat,lon,w,pi
     .      ,rottmp(3,3)

      pi=4.d0*datan(1.d0)

      if(dirn.eq.1)then
c  cartesian to spherical
        ev_sph(3) = dsqrt(ev_xyz(1)**2+ev_xyz(2)**2+ev_xyz(3)**2)
        ev_sph(1) = dasin(ev_xyz(3)/ev_sph(3))
        ev_sph(2) = datan2(ev_xyz(2),ev_xyz(1))
        lat = ev_sph(1)
        lon = ev_sph(2)
        w   = ev_sph(3)
      else
c  spherical to cartesian   
        lat = ev_sph(1)
        lon = ev_sph(2)
        w   = ev_sph(3)
        ev_xyz(1) = w*dcos(lat)*dcos(lon) 
        ev_xyz(2) = w*dcos(lat)*dsin(lon) 
        ev_xyz(3) = w*dsin(lat)            
      endif

c      print*,' in carsph lat long w ',lat*180.d0/pi
c     .         ,lon*180.d0/pi,w*180.d0/pi
c  now the jacobian matrix from spherical to cartesian
        rotmat1(1,1) = -w*dsin(lat)*dcos(lon)
        rotmat1(1,2) = -w*dcos(lat)*dsin(lon)
        rotmat1(1,3) = dcos(lat)*dcos(lon)

        rotmat1(2,1) = -w*dsin(lat)*dsin(lon)
        rotmat1(2,2) = w*dcos(lat)*dcos(lon)
        rotmat1(2,3) = dcos(lat)*dsin(lon)
  
        rotmat1(3,1) = w*dcos(lat)
        rotmat1(3,2) = 0.d0
        rotmat1(3,3) = dsin(lat)


c  now invert the jacobian to form the matrix to use from XYZ to SPH  
        do i = 1,3
          do j=1,3
            rottmp(i,j) = rotmat1(i,j) 
          enddo
        enddo
        call invert(rottmp,rotmat2,3,3)
c        call printmat(rotmat2,3,3,'rotmat2 ')
c        pause
        return
        end

      subroutine solnupd(nsit,nplates,ev_apr,tmpev_xyz,ev_xyz
     .                    ,plates,ev_sph,covsoln,covsph)

c  subroutine to add the solution corrections to the a priori estimates of the Euler vectors
c  and to convert the whole thing to spherical coordinates (including the VCV). 
c
c  P. Tregoning
c  13 October 2000

      implicit none

      integer uout,nsit,nplates,ierr,i,j,k
      real*8 ev_apr(nplates,3),tmpev_xyz(nplates*3),ev_xyz(nplates,3)
     .      ,ev_sph(nplates,3)
     .      ,covsoln(nplates*3,nplates*3),covsph(nplates*3,nplates*3) 
     .      ,lat,long,w,w2,pi,rt(3,3)
     .      ,rotmat1(3,3),tmpxyz(3),tmpsph(3),tmp1(3,3),tmpcov(3,3)
     .      ,rotmat2(3,3),smaj,smin,az    
      character plates(nplates)*20

      pi=4.d0*datan(1.d0)

c first, compute xyz a priori Euler vectors
      do i=1,nplates     
        tmpsph(1) = ev_apr(i,1)*pi/180.d0
        tmpsph(2) = ev_apr(i,2)*pi/180.d0
        tmpsph(3) = ev_apr(i,3)*pi/180.d0   
        call carsph(tmpxyz,tmpsph,2,rotmat1,rotmat2)

        ev_xyz(i,1) = tmpxyz(1)
        ev_xyz(i,2) = tmpxyz(2)
        ev_xyz(i,3) = tmpxyz(3)

c  convert the a priori Euler vectors into xyz rotations and add the corrections
        ev_xyz(i,1) = ev_xyz(i,1)  + tmpev_xyz(i*3-3+1)
        ev_xyz(i,2) = ev_xyz(i,2)  + tmpev_xyz(i*3-3+2)
        ev_xyz(i,3) = ev_xyz(i,3)  + tmpev_xyz(i*3-3+3)

c  now convert the updated xyz rotations into a spherical euler vector
        do j=1,3
          tmpxyz(j) = ev_xyz(i,j)
        enddo
        call carsph(tmpxyz,tmpsph,1,rotmat1,rotmat2)
        ev_sph(i,1) = tmpsph(1)*180.d0/pi
        ev_sph(i,2) = tmpsph(2)*180.d0/pi
        ev_sph(i,3) = tmpsph(3)*180.d6/pi
c        print*,'ev_sph = ',(ev_sph(i,j),j=1,3)

c  now compute the VCV of the spherical representation. First extract the 3x3 vector
c  for this pole
        do j=1,3  
          do k=1,3
            tmpcov(j,k) = covsoln(i*3 - 3 + j,i*3 - 3 + k) 
          enddo
        enddo  

c                                                 -1
c                                  -1            t
c  multiply it out   VCV     =  jac   VCV     Jac    where Jac = jacobian SPH -> XYZ
c                       sph               xyz     
  
        call transp(rotmat2,rt,3,3)
        call matmult(tmpcov,rt,tmp1,3,3,3)
        call matmult(rotmat2,tmp1,tmpcov,3,3,3)

c        write(*,100)plates(i),(ev_sph(i,j),j=1,3)
c     .       ,dsqrt(tmpcov(1,1))*180.d0/pi
c     .       ,dsqrt(tmpcov(2,2))*180.d0/pi,dsqrt(tmpcov(3,3))*180.d6/pi
c     .       ,tmpcov(1,2)/dsqrt(tmpcov(1,1)*tmpcov(2,2))  
c
c100     format(a,3f8.3,2f6.2,f7.4,f10.6)

c  convert it to semi-maj and semi-min error ellipse
        call errell(tmpcov(1,1),tmpcov(2,2),tmpcov(1,2),smaj,smin,az)

c  store it away
        do j=1,3  
          do k=1,3
            covsph(i*3-3+j,k) = tmpcov(j,k)
          enddo
        enddo
                                 
      enddo

c      call printmat(covsph,nplates*3,3,'covsph1 ')

      return
      end

      subroutine errell(vlat,vlong,vlatlong,smaj,smin,az)

c  subroutine to compute the semi-maj, semi-min axes and the azimuth of the
c  semi-major axis. It will use the VCV of the lat, long and latlong
c
c  P. Tregoning
c  17 October 2000

      implicit none

      real*8 vlat,vlong,vlatlong,smaj,smin,theta,tmp,pi,theta2,az
      integer i,j

      pi=4.d0*datan(1.d0)


c  now convert back to semi-maj, semi-min and orientation ellipse
          theta = 0.5*datan(2*vlatlong/(vlat-vlong))
          theta2 = theta + pi/2.d0

          smaj = dsqrt(vlat*dcos(theta)**2 + 
     .                    vlong*dsin(theta)**2 + 
     .                    vlatlong*dsin(2*theta))  

          smin = dsqrt(vlat*dcos(theta2)**2 + 
     .                    vlong*dsin(theta2)**2 + 
     .                    vlatlong*dsin(2*theta2))  


c check which of the two angles give the semi-major and which the semi-minor
           if(smin.gt.smaj)then
             az = theta2
             tmp = smaj
             smaj = smin
             smin = tmp
           else
             az = theta
           endif

c  convert mathematical theta into surveying brg from north      
          az = az*180.d0/pi
          if(az.lt.0)az = az+360.d0

c          write(*,120)smaj*180.d0/pi,smin*180.d0/pi,az
c120       format(5x,'semi-major axis: ',f6.3,/
c     .          ,5x,'semi-minor axis: ',f6.3,/
c     .          ,5x,'bearing of majr: ',f7.3,//)
          
      return
      end

      subroutine calcvel(xyz,nplates,iplate,nsit,sit,lat,long,omega
     .                  ,cvelxyz)
                            
c  subroutine to calculate an a priori velocity from the a priori Euler vector
c  Harvey (1985) which states that, provided the magnitude of rotation about
c  an Euler pole is small, the rotation can be described as
c
c
c              |   1  -Nw    Mw  |
c              |                 |
c    RE(w)  =  |  Nw    1   -Lw  |
c              |                 |
c              | -Mw   Lw    1   |
c              |                 |
c              
c              
c
c
c  where N = sin(lat)
c       
c        L = cos(lat)cos(long)
c
c        M = cos(lat)sin(long)
c
c        w = angular rotation about Euler Pole at location (lat, long)
c
c  Written by Paul Tregoning
c
c  10 October 2000
c
c

      implicit none

      integer nsit,nplates,sit,iplate,i,j
      real*8 xyz(nsit,3),lat,long,omega,cvelxyz(3)
      real*8 l,m,n,pi,R(3,3),xyz1(3),xyz2(3)

      pi=4.d0*datan(1.d0)
      
      
c zero out the starting arrays
      do i=1,3
        do j=1,3
          R(i,j) = 0.d0
        enddo
        R(i,i) = 1.d0
      enddo


c  put the xyz coordinates into a 3x1 matrix 
      do i=1,3
        xyz1(i) = xyz(sit,i)
      enddo

c  calculate the direction cosines of the Euler vector
      n = dsin(lat*pi/180.d0)
      m = dcos(lat*pi/180.d0)*dsin(long*pi/180.d0)
      l = dcos(lat*pi/180.d0)*dcos(long*pi/180.d0)

c  now form up the rotation matrix
      R(1,2) = -1.d0 * n * (omega*pi/180.d0)       
      R(1,3) =  m * (omega*pi/180.d0)
      R(2,1) = -1.d0*R(1,2)
      R(2,3) = -1.d0 * l * (omega*pi/180.d0)
      R(3,1) = -1.d0 * R(1,3)
      R(3,2) = -1.d0 * R(2,3)

c then the coordinate rotation is simply this matrix * the XYZ coords
      call matmult(R,xyz1,xyz2,3,3,1)


c  now calculate the difference in the coordinates to get the linear change in position  
      do i=1,3
        cvelxyz(i) = xyz2(i) - xyz1(i)
      enddo

c      write(*,'(a,3f10.4)')' calculated velocity is',cvelxyz
      end

      subroutine printmat(mat,m,n)

c  prints out a m x n matrix
c  PT 9/6/92

      implicit none
      integer m,n,i,j
      real*8 mat(m,n)
      character*1 ent
        
      do 10 i=1,m
          write(*,*)(mat(i,j),j=1,n)
10    continue

      print*,' '
c     write(*,'(a)')' press return to continue : '
c     read(*,'(a)')ent

      return
      end
