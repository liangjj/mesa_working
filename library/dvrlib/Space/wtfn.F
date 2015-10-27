*deck wtfn.f 
c***begin prologue     wtfn
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           wtfn functions
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            weight functions and their first and second
c***                   derivatives
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       wtfn
      subroutine wtfn(px,pwt,type,alpha,beta,ngot,n)
c
      implicit integer (a-z)
      character*(*) type
      real*8 x, wt, alpha, beta
#ifdef DECPOINTER
      integer*8 px, pwt
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*8 px, pwt
#endif SGIPOINTER
      common/io/inp, iout      
      pointer (px,x(1))
      pointer (pwt,wt(1))
      need=wptoin(3*n)
      call memory(need,pwt,ngot,'weights',0)
      p2wt=1
      p2dwt=p2wt+n
      p2ddwt=p2dwt+n
      if(type.eq.'legendre') then
         call vfill(wt(p2wt),1.d0,n)
         call rzero(wt(p2dwt),n)
         call rzero(wt(p2ddwt),n)
      elseif(type.eq.'hermite') then
         call hermit(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,n)
      elseif(type.eq.'laguerre') then
         call lagure(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,alpha,n)        
      elseif(type.eq.'chebyshev-1') then
         call cheb1(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,n)        
      elseif(type.eq.'chebyshev-2') then
         call cheb2(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,n)        
      elseif(type.eq.'jacobi') then
         call cheb2(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,alpha,beta,n)        
      else
         call lnkerr('quit. error in weight type')
      endif
      return
      end
















