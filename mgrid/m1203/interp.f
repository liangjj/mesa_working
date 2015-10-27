*deck interp
c***begin prologue     interp
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            interpolate coarse grid solution to fine grid
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       interp
      subroutine interp(uf,uc,nf,nc)
c
      implicit integer (a-z)
      real*8 uc, uf
      character*80 title
      dimension uc(nc,nc,nc), uf(nf,nf,nf)
      common/io/inp, iout
c     fill common coarse/fine points      
      do 10 kc=1,nc
         kf=2*kc-1
         do 20 jc=1,nc
            jf=2*jc-1
            do 30 ic=1,nc
               if=2*ic-1
               uf(if,jf,kf) = uc(ic,jc,kc)
   30       continue
   20    continue   
   10 continue
c    do all the points in the (x,y) planes containing the original
c    coarse grid points.  for each of these planes the value of z
c    is fixed and must change by two in order to make room for the planes
c    in between which will be filled in later. 
      do 40 kf=1,nf,2
         do 50 jf=1,nf,2
c        we are interpolating along x         
            do 60 if=2,nf-1,2
               uf(if,jf,kf) = .5d0*( uf(if+1,jf,kf) + uf(if-1,jf,kf) )
   60       continue
   50    continue   
         do 70 jf=2,nf-1,2
c        we are interpolating along y         
            do 80 if=1,nf
               uf(if,jf,kf) = .5d0*( uf(if,jf+1,kf) + uf(if,jf-1,kf) )
   80       continue   
   70    continue
   40 continue
c     now fill in the (x,y) planes in between the original coarse grid points.
      do 90 if=1,nf
         do 100 jf=1,nf
            do 200 kf=2,nf-1,2
               uf(if,jf,kf) = .5d0*( uf(if,jf,kf+1) + uf(if,jf,kf-1) )
  200       continue
  100    continue
   90 continue  
      return
      end





