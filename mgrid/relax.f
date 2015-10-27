*deck relax
c***begin prologue     relax
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            red/black gauss-seidel iteration
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       relax
      subroutine relax(u,rhs,h,k,n,type)
c
      implicit integer (a-z)
      real*8 u, rhs, h, k, hsq, ksq, fac
      character*(*) type
      character*80 title
      dimension u(n,n), rhs(n,n)
      common/io/inp, iout
c      title='input rhs'
c      call prntrm(title,rhs,n,n,n,n,iout)
c      title='initial solution'
c      call prntrm(title,u,n,n,n,n,iout)  
      hsq=h*h
      ksq=k*k*hsq
      fac=1.d0/(4.d0-ksq) 
      if(type.eq.'gauss-seidel') then
         jsw=1    
         do 10 ipass=1,2
            isw=jsw
            do 20 j=2,n-1
               do 30 i=isw+1,n-1,2
                  u(i,j)=fac*( u(i+1,j) + u(i-1,j) + u(i,j+1) 
     1                                  + u(i,j-1)
     2                                  + rhs(i,j) )
 30            continue
               isw=3-isw
 20         continue
            jsw=3-jsw   
 10      continue
      elseif(type.eq.'jacobi') then
            do 40 j=2,n-1
               do 50 i=2,n-1
                  u(i,j)=fac*( u(i+1,j) + u(i-1,j) + u(i,j+1) 
     1                                  + u(i,j-1)
     2                                  + rhs(i,j) )
 50            continue
 40         continue
       else
            call lnkerr('error in relaxation definition')
       endif   
c      title='relaxing'
c      call prntrm(title,u,n,n,n,n,iout)  
      return
      end





