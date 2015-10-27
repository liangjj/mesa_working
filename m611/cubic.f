*deck @(#)cubic.f	5.2 11/28/95
      subroutine cubic(a1,a2,a3,root,nreal)
c***begin prologue     cubic.f
c***date written       yymmdd  
c***revision date      11/28/95      
c
c***keywords           
c***author             
c***source             @(#)cubic.f	5.2   11/28/95
c***purpose            returns roots of cubic polynomial 
c***description
c   given the polynomial x**3 +a1*x**2 +a2**x +a3 =0
c   return the number of real roots, and their magnitude.
c
c***references
c                      w.h. press, b.p. flannery, s.a. teukolsky, 
c                      and w.t. vetterling, "numerical recipes",
c                      cambridge university press, p. 146,1987.
c
c***routines called
c
c***end prologue       cubic.f
      implicit none
c     --- input variables -----
      real*8 a1,a2,a3
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 root(3)
c     --- output variables ---
      integer nreal
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      real*8 q,r,pi
      real*8 zero,one
      real*8 two,three,four,nine,twent7,fift4
      real*8 theta,acos,cos,sqrt,sign,abs,atan
      real*8 t1,t2,t3
c
      parameter (zero=0.0d+00,two=2.0d+00)
      parameter (three=3.0d+00,nine=9.0d+00)
      parameter (twent7=27.0d+00,fift4=54.0d+00)
      data one/1.0d0/,four/4.0d0/
c
      common/io/inp,iout
c
c     --- initialize the solutions
      pi=four*atan(one)
      call rzero(root,3) 
c
      q=(a1*a1-three*a2)/nine
      r=(two*(a1**3) -nine*a1*a2 +twent7*a3)/fift4
c
c     --- test for the number of real roots
      if((q**3-r*r).ge.zero) then
c        three real roots
         nreal=3
         theta=acos(r/sqrt(q**3))
         root(1)=-two*sqrt(q)*cos(theta/three) - a1/three
         root(2)=-two*sqrt(q)*cos((theta+two*pi)/three) - a1/three
         root(3)=-two*sqrt(q)*cos((theta+four*pi)/three) -a1/three
      else
c        one real root
         nreal=1
         t1=sqrt(r**2-q**3)+abs(r)
         t1=t1**(one/three)
         t2=q/t1
         t3=t1+t2
         t3=-sign(one,r)*t3
         root(1)=t3 -a1/three
      endif
c
c
      return
      end
