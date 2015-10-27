*deck resinj
c***begin prologue     resinj
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            calculate residual and inject to coarse grid
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       resid
      subroutine resinj(uf,f,res,uc,h,k,nf,nc)
c
      implicit integer (a-z)
      real*8 uf, f, res, uc, h, k, rms
      character*80 title
      dimension uf(nf,nf), f(nf,nf), res(nf,nf), uc(nc,nc)
      common/io/inp, iout
      call resid(uf,f,res,rms,h,k,nf,'both')
      do 10 i=2,nc-1
         i2=i+i-1
         do 20 j=2,nc-1
            j2=j+j-1
            uc(i,j)=2.d0*res(i2,j2)
 20      continue   
 10   continue
c      title='coarse solution after injection'
c      call prntrm(title,uc,nc,nc,nc,nc,iout)
      return
      end





