c $Header: scmphi.f,v 1.2 92/12/24 15:30:50 bis Exp $
*deck scmphi.f
      subroutine scmphi(cphi,sphi,phifn,n,mu)
      implicit integer (a-z)
      common /io/ inp, iout
      dimension cphi(n), sphi(n), phifn(n)
      real *8 cphi, sphi, phifn, const
      real *8 pi, twopi
      complex*16 cmp
      data pi / 3.14159265358979323846d+00 /
      twopi=2.d0*pi
      const=sqrt(1.d0/twopi)
      if(mu.eq.0) then
         do 10 i=1,n
            phifn(i)=const
   10    continue
      else
         const=const*sqrt(2.d0)
         absmu=abs(mu)
         if (mu.gt.0) then
             do 20 i=1,n
                cmp=(dcmplx(cphi(i),sphi(i)))**absmu
                phifn(i)=real(cmp)*const
   20        continue
         else
             do 30 i=1,n
                cmp=(dcmplx(cphi(i),sphi(i)))**absmu                         
                phifn(i)=imag(cmp)*const
   30        continue
         endif
      endif
      return
      end
