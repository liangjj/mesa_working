c $Header: miscfn.f,v 1.1 92/12/23 17:37:45 bis Exp $
*deck miscfn.f
      subroutine miscfn(phi,cphi,sphi,nphi)
      implicit integer (a-z)
      real *8 phi, cphi, sphi
      dimension phi(nphi), cphi(nphi), sphi(nphi)
      do 10 pnt=1,nphi
         cphi(pnt)=cos(phi(pnt))
         sphi(pnt)=sin(phi(pnt))
   10 continue
      return
      end
