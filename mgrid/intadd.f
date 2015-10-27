*deck intadd
c***begin prologue     intadd
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            interpolate to fine grid
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       intadd
      subroutine intadd(uc,uf,nc,nf)
c
      implicit integer (a-z)
      real*8 uf, uc
      character*80 title
      dimension uf(nf,nf), uc(nc,nc)
      common/io/inp, iout      
      do 10 i=2,nc
         i2=i+i-1
         do 20 j=2,nc
            j2=j+j-1
            uf(i2,j2-1) = .5d0*( uc(i,j) + uc(i,j-1) ) + uf(i2,j2-1)
            uf(i2-1,j2) = .5d0*( uc(i,j) + uc(i-1,j) ) + uf(i2-1,j2)
 20      continue
 10   continue    
      call rzero(uc,nc*nc)
c      title='fine grid solution by interpolation'
c      call prntrm(title,uf,nf,nf,nf,nf,iout)
      return
      end





