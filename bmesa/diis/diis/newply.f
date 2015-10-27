*deck newply.f
c***begin prologue     newply
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            evaluate dvr polynomials at dvr points.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       newply
      subroutine newply(p,dp,ddp,a,b,pn,dpn,ddpn,q,u,left,right,
     1                  nleft,nright,n)
      implicit integer (a-z)
      real*8 p, dp, ddp, a, b, pn, dpn, ddpn, q, u
      real*8 left, right
      dimension p(n,0:n-1), dp(n,0:n-1), ddp(n,0:n-1)
      dimension a(0:n), b(0:n) 
      dimension pn(n,0:n-1), dpn(n,0:n-1), ddpn(n,0:n-1)
      dimension q(n), u(n,n)
      common/io/inp, iout 
      call gpoly(p,dp,ddp,q,a,b,left,right,nleft,nright,n,n,.false.)
c     transform the polynomials and their derivatives to the new
c     representation      
      call ebc(pn(1,0),p(1,0),u,n,n,n)
      call ebc(dpn(1,0),dp(1,0),u,n,n,n)
      call ebc(ddpn(1,0),ddp(1,0),u,n,n,n)
      return
      end       
