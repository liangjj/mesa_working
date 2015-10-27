*deck cmpare
c***begin prologue     cmpare
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            compare exact and approximate solution and calculate
c***references         various residuals
c
c***routines called    iosys, util and mdutil
c***end prologue       cmpare
      subroutine cmpare(u,ex,f,h,k,n)
c
      implicit integer (a-z)
      real*8 u, ex, f, h
      real*8 emax, esum, rmax, rsum, er, r, eh, rh, k, hsq, hsqi, fac
      character*80 title
      dimension u(n,n), ex(n,n), f(n,n)
      common/io/inp, iout
      hsq=h*h
      hsqi=1.d0/hsq
      fac=4.d0-k*k*hsq
      emax=0.d0
      esum=0.d0
      rmax=0.d0
      rsum=0.d0
      do 10 i=2,n-1
         do 20 j=2,n-1
            er=abs ( u(i,j)- ex(i,j) )
            if (er.gt.emax) then
                emax=er
            endif
            r=abs ( f(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + 
     1              u(i,j-1) - fac*u(i,j) )
            if (r.gt.rmax) then
                rmax=r
            endif
            esum=esum+er*er
            rsum=rsum+r*r
 20      continue    
 10   continue   
      eh=sqrt ( h*esum )
      rh=sqrt ( h*rsum )
      write(iout,1) eh, rh
      write(iout,2) emax, rmax
      call vsub(ex,ex,u,n*n)
c      title='difference between exact and numerical solution'
c      call prntrm(title,ex,n,n,n,n,iout)
 1    format(/,1x,'step size weighted rms error of solution = ',e15.8,
     1       /,1x,'step size weighted rms error of residual = ',e15.8)
 2    format(/,1x,'largest absolute difference = ',e15.8,/,1x,
     1            'largest absolute residual   = ',e15.8)
      return
      end





