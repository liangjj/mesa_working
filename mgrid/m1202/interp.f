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
      dimension uc(nc,nc), uf(nf,nf)
      common/io/inp, iout
      do 10 jc=1,nc
         jf=2*jc-1
         do 20 ic=1,nc
            uf(2*ic-1,jf) = uc(ic,jc)
 20      continue   
 10   continue
      do 30 jf=1,nf,2
         do 40 if=2,nf-1,2
            uf(if,jf) = .5d0*( uf(if+1,jf) + uf(if-1,jf) )
 40      continue   
 30   continue
      do 50 jf=2,nf-1,2
         do 60 if=1,nf
            uf(if,jf) = .5d0*( uf(if,jf+1) + uf(if,jf-1) )
 60      continue   
 50   continue
      return
      end





