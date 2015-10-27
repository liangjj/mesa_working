*deck scmphi.f
c***begin prologue     scmphi
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scmphi, link 6201, legendre functions
c***author             schneider, barry (lanl)
c***source             m6201
c***purpose            calculate sin(m*phi) and cos(m*phi)
c***references         none
c
c***routines called
c***end prologue       scmphi
      subroutine scmphi(phi,phifn,n,m)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 phi, phifn, const
      real *8 pi, twopi
      dimension phifn(n,2), phi(n)
      data pi / 3.14159265358979323846d+00 /
      twopi=2.d0*pi
      const=sqrt(1.d0/twopi)
      if(m.eq.0) then
         do 10 i=1,n
            phifn(i,1)=const
   10    continue
      else
         const=const*sqrt(2.d0)
         do 20 i=1,n
            phifn(i,1)=const*cos(m*phi(i))
            phifn(i,2)=const*sin(m*phi(i))  
   20    continue
      endif
      return
      end
