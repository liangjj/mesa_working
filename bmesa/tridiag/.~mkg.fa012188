*deck mkg.f
c***begin prologue     mkg
c***date written       950720 (yymmdd)
c***revision date             (yymmdd)
c***keywords           
c***author             schneider, barry(nsf)
c***source             @(#)util
c***purpose            set up g function for numerov
c***                   integration of one-dimensional
c***                   schroedinger equation.
c***                                         ''
c***                   the equation here is y  + f(x) y = g
c***                   and we are computing g.
c***routines called  
c***end prologue
      subroutine mkg(g,v,r,energy,l,n,zro,prnt)
      implicit integer (a-z)
      real*8  g, v, r, j, jp, jpp, y, yp, ypp
      real*8 energy, k, rr
      logical prnt, zro
      character*80 title
      dimension g(n), v(n), r(n)
      common /io/ inp, iout
      call rzero(g,n)
      if(.not.zro) then
          k=sqrt(energy)
          do 10 i=1,n
             rr=r(i)*k
             call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)
             g(i)=2.d0*v(i)*j
   10     continue
      endif      
      if (prnt) then
          title='inhomogeneity on grid'
          call prntrm(title,g,n,1,n,1,iout)
      endif
      return
      end



