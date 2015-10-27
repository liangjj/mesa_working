*deck @(#)scmphi.f	1.2  10/27/94
c***begin prologue     scmphi
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scmphi, link 6203, legendre functions
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            calculate sin(m*phi) and cos(m*phi)
c***references         none
c
c***routines called
c***end prologue       scmphi
      subroutine scmphi(cphi,sphi,phifn,n,m)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 cphi, sphi, phifn, const
      real *8 pi, twopi
      complex*16 cmp
      dimension cphi(n), sphi(n), phifn(n,2)
      data pi / 3.14159265358979323846d+00 /
      twopi=2.d0*pi
      const=sqrt(1.d0/twopi)
      if(m.eq.0) then
         do 10 i=1,n
            phifn(i,1)=const
   10    continue
      else
c        the functions for m.ne.0 have different norm (1/sqrt(pi)).
         const=const*sqrt(2.0d0)
         do 20 i=1,n
            cmp=(dcmplx(cphi(i),sphi(i)))**m
            phifn(i,1)=real(cmp)*const
            phifn(i,2)=imag(cmp)*const
   20    continue
      endif
      return
      end
