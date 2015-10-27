*deck addint
c***begin prologue     addint
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            interpolate to fine grid and add to fine grid
c***                   solution
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       addint
      subroutine addint(uf,uc,res,nf,nc)
c
      implicit integer (a-z)
      real*8 uf, uc, res
      character*80 title
      dimension uf(nf,nf), uc(nc,nc), res(nf,nf)
      common/io/inp, iout      
      call interp(res,uc,nf,nc)
      do 10 j=1,nf
         do 20 i=1,nf
            uf(i,j) = uf(i,j) + res(i,j)
 20      continue
 10   continue    
      return
      end





