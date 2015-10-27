*deck hyprad.f
c***begin prologue     hyprad
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            hyperspherical radial functions.
c***                   
c***                   R(rho,m) = J(rho,2*m+2) / (rho**2)
c***                   I(rho,m) = N(rho,2*m+2) / (rho**2)
c***
c***                   J and N are Bessel and Neumann functions.
c***
c***                   J(rho,m) goes like 
c***                   sqrt(2./{pi*rho}) * cos ( rho - pi/2 * { m + 1/2 } )
c***                   for large rho.  The other goes sine-like.
c***                   this can be written in a more standard form, with
c***                   an effective angular momentum of l = 2*m + 3/2
c***                   by using the addition theorem for cosines.  then
c***                   we see that
c***                   J(rho,m) goes like
c***                   sqrt(2./{pi*rho}) * sin ( rho - l*pi/2 )
c***
c***                   similarly, N(kr,m) can be shown to be,
c***                   - sqrt(2./{pi*rho}) * cos ( rho - l*pi/2 )
c***
c***                   thus the entire radial wavefunction has a
c***                   1/r**(5/2) type inverse radial factor.
c***references         
c
c***routines called    
c***end prologue       hyprad
      subroutine hyprad(r,jmk,djmk,nmk,dnmk,k,jb,djb,yb,dyb,fact,
     1                  m,n1,n2,top,type)
      implicit integer (a-z)
      real*8 r, jmk, djmk, nmk, dnmk 
      real*8 k, jb, djb, yb, dyb, rho, fact
      character*(*) type
      dimension r(n2,n1), jmk(n2,n1), nmk(n2,n1)
      dimension djmk(n2,n1), dnmk(n2,n1)
      dimension jb(0:*), djb(0:*), yb(0:*), dyb(0:*), fact(0:*)
      common/io/inp, iout
c
      do 10 i=1,n1
         do 20 j=1,n2
            rho=k*r(j,i)
            call radial(rho,jmk(j,i),djmk(j,i),nmk(j,i),dnmk(j,i),
     1                  jb,djb,yb,dyb,fact,top,m,type)
 20      continue   
 10   continue   
c
c    make derivative wrt r not rho
c
      call vscale(djmk,djmk,k,n1*n2)
      call vscale(dnmk,dnmk,k,n1*n2)
      return
      end       















