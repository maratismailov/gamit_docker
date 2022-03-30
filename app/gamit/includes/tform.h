c Set dimensions and store until numbers for TFORM 

c    Station and coordinate dimensions for local arrays

      integer nsdim,ndim              

      parameter (nsdim = 3000)

      parameter (ndim = nsdim*3 )
                      
c    Unit numbers, set by tform or convgeo

      integer iterm,iscrn,iprnt,idatum 
      common /tfrmcom/ iterm,iscrn,iprnt,idatum

