c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Chpsi}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck chpsi.f
***begin prologue     chpsi
***date written       011017   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development
***author             schneider, barry (nsf)
***source             
***purpose            right hand side of inhomogeneous time-dependent
***                   hamiltonian.
***description        matrix element of negative of ( i d/dt - H ) times
***                   initial state.  
***                   n = nspace = nx one dimension
***                              = ny*nx two dimensions
***                              = nz*ny*nx three dimensions
***                   note that the vectors are stored exactly as above.
***references         
***routines called    
***end prologue       chpsi
      subroutine chpsi(h1,h2,h3,driver,psi0,v,typke,n,nd,ndiff,dim,prn)
      implicit integer (a-z)
      real*8 h1, h2, h3, v
      complex*16 psi0, driver, cdum
      character*(*) typke
      logical prn
      character*80 title
      dimension nd(dim), ndiff(dim)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2)), h3(nd(3),nd(3))
      dimension driver(n), v(n), psi0(n)
      common/io/inp, iout

c     Zero the driver array

      call czero(driver,n)

c     Multiply by the diagonal potential

      call cvmul(driver,psi0,v,cdum,n,'real')

c     Now do the one body operators.

      if(dim.eq.1) then
         kb1 = ndiff(1)/2 + 1
         call h1v(h1,h1,driver,psi0,nd(1),kb1,1,typke)
      elseif(dim.eq.2) then
         kb1 = ndiff(1)/2 + 1
         kb2 = ndiff(2)/2 + 1
         call h2v(h1,h2,h1,h2,driver,psi0,nd(1),nd(2),
     1            kb1,kb2,1,typke)
      elseif(dim.eq.3) then
         kb1 = ndiff(1)/2 + 1
         kb2 = ndiff(2)/2 + 1
         kb3 = ndiff(3)/2 + 1
         call h3v(h1,h2,h3,h1,h2,h3,driver,psi0,nd(1),nd(2),nd(3),
     1            kb1,kb2,kb3,1,typke)
      endif   
      if(prn) then
         title='driving term'
         call prntcm(title,driver,n,1,n,1,iout)
      endif               
      return
 1    format(/,1x,'test sum = ',e15.8,1x,'i = ',i3,1x,'j = ',i3)
      end
