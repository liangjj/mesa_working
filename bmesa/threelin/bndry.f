*deck bndry.f
c***begin prologue     bndry
c***date written       921206   (yymmdd)
c***revision date               (yymmdd)
c***keywords           numerov
c***author             schneider, barry(lanl)
c***source             onelin
c***purpose            set boundary condition for solution of finite
c***                   difference equation.
c***description        the equation here is y  + f(x) y = g
c***
c***                   first derivative formula needed for 
c***                   boundary condition is:
c***                   ( -3*y(0) + 4*y(1) -y(2) )/(2*stp)
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis. it is well
c***                   described in Kopal's book on numerical analysis.
c***routines called  
c***end prologue
      subroutine bndry(diag,sudiag,spdiag,f,g,rhs,x,s4,stp,energy,m,
     1                  last,nfd,bcond,value)
      implicit integer (a-z)
      real *8  diag, sudiag, spdiag, f, g, rhs, stp, value, energy, k
      real *8 x, fac, fac1, s4
      character*(*) bcond 
      dimension diag(0:m), sudiag(0:m), spdiag(0:m), f(0:m), g(0:m)
      dimension rhs(0:m), x(0:m)
      common /io/ inp, iout
      do 10 i=1,m-1
         rhs(i) = stp*stp*( g(i+1) +10.d0*g(i) +g(i-1) )
   10 continue
      nfd=last-1
      s4=rhs(nfd)
      if (bcond.eq.'function') then
          rhs(nfd)=rhs(nfd)-spdiag(nfd)*value
          return 
      elseif (bcond.eq.'derivative') then
          sudiag(nfd)=sudiag(nfd)-spdiag(nfd)/3.d0
          diag(nfd)=diag(nfd)+4.d0*spdiag(nfd)/3.d0
          return 
      elseif (bcond.eq.'log-derivative') then
          k=sqrt(energy)
          value=-k*tan(k*x(last))
          fac=1.d0-stp*value-f(last)*stp*stp*.5d0
          fac1=stp*stp*g(last)*.5d0
          diag(nfd)=diag(nfd)+spdiag(nfd)/fac
          rhs(nfd)=rhs(nfd)+spdiag(nfd)*fac1/fac          
          return 
      elseif (bcond.eq.'none') then
          return 
      endif 
      return
      end



