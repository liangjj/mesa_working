c $Header: scmphi.f,v 1.2 92/12/24 15:30:50 bis Exp $
*deck scmphi.f
      subroutine scmphi(cphi,sphi,phifn,n,mu)
      implicit integer (a-z)
      common /io/ inp, iout
      dimension cphi(n), sphi(n), phifn(n,2)
      real *8 cphi, sphi, phifn, const
      real *8 pi, twopi
      complex*16 cmp
      data pi / 3.14159265358979323846d+00 /
      twopi=2.d0*pi
      const=sqrt(1.d0/twopi)
      if(mu.eq.0) then
         do 10 i=1,n
            phifn(i,1)=const
   10    continue
         return
      else
         const=const*sqrt(2.d0)
         do 20 i=1,n
            cmp=(dcmplx(cphi(i),sphi(i)))**mu
            phifn(i,1)=real(cmp)*const
            phifn(i,2)=imag(cmp)*const
   20    continue
         return
      endif
      return
      end
